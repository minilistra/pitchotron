/*  Copyright (C) 2007 Karel 'Clock' Kulhavy - Twibright Labs
 *  This is a part of the Ronja project. pix-y4m a tool to convert
 *  BRL-CAD .pix files into the YUV4MPEG2 format (a 4:2:0 Y'CbCr with
 *  extra simple ASCII headers).
 *  
 *  This is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This software is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this software; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
 *  USA.
 */

/* http://www.poynton.com/PDFs/ColorFAQ.pdf */

#define _GNU_SOURCE
/* According to
 * http://www2.informatik.uni-halle.de/lehre/c/libc/libc.texinfo_1.html#SEC12
 * "_GNU_SOURCE If you define this macro, everything is included: ANSI C,
 * POSIX.1, POSIX.2, BSD, SVID, and GNU extensions. In the cases where POSIX.1
 * conflicts with BSD, the POSIX definitions take precedence. "
 */

#include <stdlib.h> /* exit */
#include <stdio.h> /* sscanf */
#include <math.h> /* floor */
#include <string.h> /* memcpy */
#include <stdlib.h> /* exit, random */
#include <assert.h>
#include <sys/time.h> /* Gettimeofday */
#include <unistd.h> /* sleep */
#include <limits.h>
#include <errno.h>

#include "constants.h"

#define NO_SUBSAMPLE 0
#define TEST_STRIP_HALF 0

#define DEBUG_REV 0
#define DEBUG_Yc 217
#define DEBUG_Cb 173
#define DEBUG_Cr 131

#define ALGO_ORDINARY 0
#define ALGO_HYPERLUMA_1 1
#define ALGO_HYPERLUMA_2 2
#define ALGO_LUMINAPLEX 3

#define ALGO ALGO_LUMINAPLEX

#define ITIMER_USED ITIMER_VIRTUAL

unsigned input_w, input_h, output_w, output_h, rate; /* Width, height (pixels), frame rate (fps) */
unsigned short *d1; /* Data buffer. The memory organization: each pixel has
		       RGB successively. Top line first. Width, height:
		       input_w, input_h */
unsigned short *d2; /* Data buffer. The memory organization: each pixel has
		       RGB successively. Top line first. Width, height:
		       output_w, output_h */
unsigned char *d3; /*  Organization by planes - first R' plane, then G' plane,
		       then B' plane. Later holds Y'CbCr. Width, height:
		       output_w, output_h */
unsigned char *line_buffer; /* Input line buffer. Allocated to 3*w. */
unsigned char *progname=(unsigned char *)"pix-y4m";
unsigned short forward_table[256]; /* Converts from the input to photometrically
				      linear. */
unsigned char backward_table[65536]; /* Converts from photometrically linear
					to the TV standard (gamma=1/TV_GAMMA) */
unsigned short backward_table_16[65536];
double input_gamma=1.0/TV_GAMMA; /* Gamma at which the input was written.
                             Default - the TV standard */
int ppm; /* If this is set, generate PPM instead of YUV4MPEG */
unsigned temporal_resample; /* 1 means normal, 2 means average each 2
			       consecutive frames and then dump them at
			       once (halving the frame rate), etc. */
unsigned char desaturate; /* 1=desaturate the testcard */

/* Question: Give Cb (16-240) and Cr (16-240) and gamma-luminance (0-255),
 * what is the luma value to achieve this gamma-luminance with these Cb and
 * Cr? The whole table size is 225*225*255=12960000 (12.96MB) */
unsigned char hyperluma_table[225*225*256];
#ifndef NOTABLE
extern unsigned char hyperluma_table_compressed[];
extern unsigned char *hyperluma_table_compressed_end;
#endif
int hyperluma_enable=1; /* You can turn it off permanently here. */
int demo; /* Hyperluma demo - hyperluma (and the logo) are being turned
	     on and off. */
float hyperluma_switch_period=1; /* For demo purpose */
unsigned long long total_errors[5]; /* Yc, R', G', B' errors and reused
				     error*/
unsigned long long total_quadpixels; /* Total number of 2x2 pixel blocks */
unsigned long long total_iterations; /* Counts how many times the
					perceived_error function has been
					called. */
static unsigned long repetitions; /* We assume the frame doesn't have 16
				     Gipixels or more */
static unsigned long non_repetitions;
static unsigned long last_error; /* This is for reuse, to know what error
				    the last value had. */

#define SRAL(x,shift) ((((x)+0x40000000)>>(shift))-(0x40000000>>(shift)))
#define GL_FORWARD(foo) \
	if (foo<256){\
		foo=forward_table[foo];\
	}else foo=65535;
/* Y 0-255, Cb and Cr 16-240
 * Returns luminance powered to TV_GAMMA, range 0-255
 */
long gamma_luminance(int y, int cb, int cr
#if DEBUG_REV
		, int debug
#endif
		)
{
	long r,g,b;
	unsigned long luminance;

	y-=16;
	cb-=128;
	cr-=128;

	r=SRAL(LUMA_16*y+R_CR_16*cr+0x8000,16);
	g=SRAL(LUMA_16*y+G_CB_16*cb+G_CR_16*cr+0x8000,16);
	b=SRAL(LUMA_16*y+B_CB_16*cb+0x8000,16);

#if DEBUG_REV
	if (debug)
#endif
	if (r>=256||g>=256||b>=256) return 131072; /* Out of range */
	if (r<0) r=0;
	if (g<0) g=0;
	if (b<0) b=0;

	GL_FORWARD(r);
	GL_FORWARD(g);
	GL_FORWARD(b);

	luminance=(LUMI_R_16*r+LUMI_G_16*g+LUMI_B_16*b+0x8000)>>16;

	return backward_table_16[luminance];
}

/* gl = gamma-luminance */
unsigned char reverse(unsigned char gl, unsigned char cb, unsigned
		char cr)
{
	unsigned char mask, y; /* y=luma */
	long diff, diff1;
	long gl_16=((unsigned)gl<<8)+gl;
	int diff_corresponds_to=-1;
#if DEBUG_REV
	int debug;

	if (gl==DEBUG_Yc &&cb==DEBUG_Cb &&cr==DEBUG_Cr) debug=1; else debug=0;
	if (debug) fprintf(stderr,"gl_16=%ld\n", gl_16);
#endif

	for (mask=0x80,y=0;mask; mask>>=1) {
		y|=mask;
		diff=gamma_luminance(y, cb, cr
#if DEBUG_REV
				, debug
#endif
				);
#if DEBUG_REV
		if (debug) fprintf (stderr,"Y'=%02x -> gl_16=%ld\n",y, diff);
#endif
		diff-=gl_16;
		diff_corresponds_to=y;
		if (!diff) {
			return y; /* Found exactly */
		}
		if (diff>0) y-=mask;
		/* Reset the bit */
	}
	if (y>=255) return y; /* Can't test one higher */
	if (diff_corresponds_to!=y)
		diff=gamma_luminance(y, cb, cr
#if DEBUG_REV
				, debug
#endif
				)-gl_16;
	diff1=gamma_luminance(y+1, cb, cr
#if DEBUG_REV
			, debug
#endif
			)-gl_16;
	if (labs(diff)>labs(diff1)){
		/* diff1 is better */
		return y+1;
	}else{
		/* diff is better */
		return y;
	}
}

/* Returns forbidden_0_to_7 for next char.
 * Prints a char for a C string. */
int putchar_c(int c)
{
	static int char_pos=0; /* To not makes lines excessively long. */
	static int forbidden_0_to_7;
	FILE *output=stdout;

	if (char_pos>=70){
		fputs("\"\n\"",output);
		char_pos=0;
	}
	switch(c){
		case '\n':
			fputs("\\n",output);
two:
			char_pos+=2;
			forbidden_0_to_7=0;
			break;

		case '\t':
			fputs("\\t",output); goto two;

		case '\b':
			fputs("\\b",output); goto two;

		case '\r':
			fputs("\\r",output); goto two;

		case '\f':
			fputs("\\f",output); goto two;

		case '\\':
			fputs("\\\\",output); goto two;

		case '\'':
			fputs("\\\'",output); goto two;

		default:
			if (c<' '||c=='"'||c=='?'||c>126
				||(c>='0'&&c<='7'&&forbidden_0_to_7)){
				fprintf(output,"\\%o",c);
				if (c>=0100) char_pos+=4;
				else if (c>=010) char_pos+=3;
				else char_pos+=2;
				forbidden_0_to_7=(c<0100);
			}else{
				fprintf(output,"%c",c);
				char_pos++;
				forbidden_0_to_7=0;
			
			}
			break;
	}
	return forbidden_0_to_7;
}


/* bit 7: 0: lowest 2 bits value, bits 2-6 repetitions
 *        1: bits 0-6 repetitions, following byte value
 */
void write_rle(unsigned char c, unsigned long count)
{
	unsigned long eat;

	while(count){
		eat=count;
		if (!(c&~3)){
			if (eat>=31) eat=31;
			putchar_c(c|(eat<<2));
		}else{
			if (eat>=127) eat=127;
			putchar_c(eat|0x80);
			putchar_c(c);
		}
		count-=eat;
	}
}

/* bit 7: 0: lowest 2 bits value, bits 2-6 repetitions
 *        1: bits 0-6 repetitions, following byte value
 * result=-1: called at the end to flush the last RLE
 */
void hyperluma_table_rle(int result)
{
	static int last_result=0;
	static int last_derivative=0;
	static unsigned long rle_count=0;
	int derivative;

	derivative=(result-last_result)&255;
	last_result=result;
	if (derivative==last_derivative&&result!=-1){
		rle_count++;
	}else {
		write_rle(last_derivative, rle_count);
		last_derivative=derivative;
		rle_count=1;
	}

}

int calc_table(int compressed)
{
	unsigned gl, cb, cr, result;

	if (compressed) printf("const unsigned char "
		"hyperluma_table_compressed[]=\n\"");
	for (cr=16;cr<=240;cr++)
		for (cb=16;cb<=240;cb++)
			for (gl=0;gl<256;gl++){

				result=reverse(gl, cb, cr);
				if (result<16) result=16; /* Avoid sync */
				if (result>235) result=235; /* Avoid sync */

				if (compressed) hyperluma_table_rle(result);
				else putchar(result);
			}
	if (compressed){
		hyperluma_table_rle(-1);
		putchar_c(0); /* Terminator - 0 repetitions */
		printf("\";\n"
			"\nconst unsigned char "
			"*hyperluma_table_compressed_end"
			"=hyperluma_table_compressed"
			"+sizeof hyperluma_table_compressed;\n");
	}
	return 0;
}

void my_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	if (fwrite(ptr, size, nmemb, stream)
		!=nmemb){
		fprintf(stderr,"%s: Error occured writing %lu bytes"
			"into a file stream: ", progname,
			(unsigned long)size*nmemb);
		perror(NULL);
		exit(1);
	}
}

/* Converts 0-255 to BT.601 Y'CbCr */
void rgb2ycbcr(unsigned char *d, unsigned plane)
{
	unsigned char *r, *g, *b;
	int y, cb,cr;

	r=d;
	g=r+plane;
	b=g+plane;
	for (;plane;plane--){

		y=16+((LUMA_R_16**r+LUMA_G_16**g+LUMA_B_16**b+0x8000)>>16);
		cb=0x80+((CB_R_16**r+CB_G_16**g+CB_B_16**b+0x8000)>>16);
		cr=0x80+((CR_R_16**r+CR_G_16**g+CR_B_16**b+0x8000)>>16);

		*r++=y;
		*g++=cb; /* Cb */
		*b++=cr; /* Cr */
	}
}


/* Converts 0-255 r'g'b' into gamma-luminance, cb, cr */
void rgb2glcbcr(unsigned char *d, unsigned plane)
{
	unsigned char *r, *g, *b;
	unsigned linr, ling, linb;
	unsigned luminance;
	int cb,cr;

	r=d;
	g=r+plane;
	b=g+plane;
	for (;plane;plane--){
		linr=forward_table[*r];
		ling=forward_table[*g];
		linb=forward_table[*b];

		luminance=(LUMI_R_16*linr+LUMI_G_16*ling
				+LUMI_B_16*linb+0x8000)>>16;
		/* Checked against overflow. Luminance is now 0-65535. */

		cb=0x80+((CB_R_16**r+CB_G_16**g+CB_B_16**b+0x8000)>>16);
		cr=0x80+((CR_R_16**r+CR_G_16**g+CR_B_16**b+0x8000)>>16);

		*r++=backward_table[luminance]; /* Store the gamma-corrected
						   luminance */
		
		*g++=cb; /* Cb */
		*b++=cr; /* Cr */
	}
}

unsigned char hyperluma(unsigned char gl, unsigned char cb, unsigned char cr)
{
	unsigned index;
	if (cb<16) cb=16;
	if (cr<16) cr=16;
	if (cb>240) cb=240;
	if (cr>240) cr=240;
	cb-=16;
	cr-=16;
	index=gl+((cr*225+cb)<<8);
	return hyperluma_table[index];
}

void subsample_horizontal_cositing(unsigned char *ptr, unsigned len)
{
	unsigned char *endptr=ptr+len;

	*ptr=(3*ptr[0]+ptr[1]+2)>>2;
	ptr+=2; /* Avoid segfault */

	for (;ptr+1<endptr;ptr+=2){
		*ptr=(ptr[-1]+(ptr[0]<<1)+ptr[1]+2)>>2;
	}
	if (ptr<endptr) *ptr=(ptr[-1]+3*ptr[0]+2)>>2;
}

/* In each pair, averages the pair and puts the result into the left member.
 * The right member to be "interpolated" (=copied) later. */
void subsample_horizontal_interstitial(unsigned char *ptr, unsigned len)
{
	unsigned char *endptr=ptr+len;

	for (;ptr+1<endptr;ptr+=2){
		*ptr=(ptr[0]+ptr[1]+1)>>1;
	}
}

void interpolate_horizontal_cositing(unsigned char *ptr, unsigned len)
{
	unsigned char *endptr=ptr+len;

	for (ptr+=1;ptr+1<endptr;ptr+=2){
		*ptr=(ptr[-1]+ptr[1]+1)>>1;
	}
	if (ptr<endptr) *ptr=ptr[-1];
}

void interpolate_horizontal_interstitial(unsigned char *ptr, unsigned len)
{
	unsigned char *endptr=ptr+len;

	for (ptr+=1;ptr<endptr;ptr+=2){
		*ptr=ptr[-1];
	}
}

void subsample_vertical_cositing(unsigned char *ptr,
		unsigned w, unsigned h)
{
	unsigned y;
	unsigned char *endptr;

	endptr=ptr+w;

	for (endptr=ptr+w;ptr<endptr; ptr++)
		ptr[0]=(3*ptr[0]+ptr[1]+2)>>2;
	ptr+=w;
	for (y=2;y+1<h;y+=2){
		for (endptr=ptr+w;ptr<endptr; ptr++)
			ptr[0]=((ptr[0]<<1)+ptr[w]+ptr[-w]+2)>>2;
		ptr+=w;
	}
	if (y<h){
		for (endptr=ptr+w;ptr<endptr; ptr++)
			ptr[0]=(3*ptr[0]+ptr[-w]+2)>>2;
	}
}

void subsample_vertical_interstitial(unsigned char *ptr,
		unsigned w, unsigned h)
{
	unsigned y;
	unsigned char *endptr;

	for (y=0;y+1<h;y+=2){
		for (endptr=ptr+w;ptr<endptr; ptr++)
			ptr[0]=(ptr[0]+ptr[w]+1)>>1;
		ptr+=w;
	}
}

void interpolate_vertical_cositing(unsigned char *ptr, unsigned w, unsigned h)
{
	unsigned y;
	unsigned char *endptr;

	ptr+=w;
	for (y=1;y+1<h;y+=2){
		for (endptr=ptr+w;ptr<endptr;ptr++)
			ptr[0]=(ptr[-w]+ptr[w]+1)>>1;
		ptr+=w;
	}
	if (y<h) memcpy(ptr, ptr-w, w);
}

void interpolate_vertical_interstitial(unsigned char *ptr, unsigned w, unsigned h)
{
	unsigned y;

	for (y=0;y+1<h;y+=2){
		memcpy(ptr+w, ptr, w);
		ptr+=(w<<1);
	}
}

void write_frame_header(void)
{
	unsigned char frame_header[]="FRAME\n";
	my_fwrite(frame_header, sizeof(frame_header)-1, 1, stdout);
}

/* Dumps the data from d, output_w, output_h */
void dump_y4m_frame(unsigned char *d, unsigned w, unsigned h)
{
	unsigned x, y; /* Counters */
	unsigned char *rptr, *gptr, *bptr; /* R'/Cr, green/get, B'/Cb pointer */
	unsigned char *wptr; /* Write ptr for moving the chroma samples
				together */
	unsigned offset=w*h;

	/* Convert the R'G'B' buffer into an Y'CbCr buffer */
	if (hyperluma_enable) rgb2glcbcr(d, offset);
	else rgb2ycbcr(d, offset);

	for (rptr=d+offset; rptr<d+3*offset; rptr+=w){

#if !NO_SUBSAMPLE
		subsample_horizontal_interstitial(rptr, w);
#endif

		if (hyperluma_enable)
			interpolate_horizontal_interstitial(rptr, w);
	}

#if !NO_SUBSAMPLE
	/* Seems according to what they write as is default for
	 * YUV4MPEG2 and how mplayer interpolates it into png output,
	 * the YUV4MPEG2 has interstitial position in both directions. */
	subsample_vertical_interstitial(d+offset, w, h);
	subsample_vertical_interstitial(d+(offset<<1), w, h);
#endif

	if (hyperluma_enable){
		interpolate_vertical_interstitial(d+offset, w, h);
		interpolate_vertical_interstitial(d+(offset<<1), w, h);
	}
	
#if TEST_STRIP_HALF
	/* Strip half of the screen of chroma */
	memset(d+offset, 128, offset>>1);
	memset(d+(offset<<1), 128, offset>>1);
#endif

	if (hyperluma_enable){
		/* Apply the hyperluma */
		for(rptr=d; rptr<d+offset; rptr++)
			*rptr=hyperluma(*rptr, rptr[offset], rptr[offset<<1]);
	}

	/* Shuffle the chroma so the data can be written out at once. */
	bptr=d+offset; /* Cb */
	rptr=bptr+offset; /* Cr */
	wptr=bptr;
	for (y=0;y<h-1;y+=2){
		gptr=bptr+y*w; /* gptr now doesn't mean green pointer,
				  but get pointer. */
		for (x=0;x<w-1;x+=2, gptr+=2)
			*wptr++=*gptr; /* Only subsampling, now the
					   decimation has already been
					   done. */
	}
	for (y=0;y<h-1;y+=2){
		gptr=rptr+y*w;
		for (x=0;x<w-1;x+=2, gptr+=2)
			*wptr++=*gptr;
	}
	/* Now the data to be written begin at d and are wptr-d long. */

	/* Write out the frame */
	write_frame_header();
	my_fwrite(d, wptr-d, 1, stdout);
}

/* Writes a line into the 'd1' buffer. Converts from R'G'B' to photometrically
 * linear RGB on the fly. The relevant dimensions are input_w and input_h. */
void write_line(unsigned line)
{
	unsigned char *iptr;
	unsigned short *optr;

	for
		(iptr=line_buffer,optr=d1+(input_h-1-line)*input_w*3;
		 	iptr<line_buffer+3*input_w;iptr+=3,optr+=3){
		*optr=forward_table[*iptr]; /* Red */
		optr[1]=forward_table[iptr[1]]; /* Green */
		optr[2]=forward_table[iptr[2]]; /* Blue */
	}
}

/* Writes a line into the 'd1' buffer. Converts from R'G'B' to photometrically
 * linear RGB on the fly. The relevant dimensions are input_w and input_h. */
void write_line_div_add(unsigned line, unsigned resample)
{
	unsigned char *iptr;
	unsigned short *optr;

	/* We do not round here because it could case arithmetic overflow. */
	for
		(iptr=line_buffer,optr=d1+(input_h-1-line)*input_w*3;
		 	iptr<line_buffer+3*input_w;iptr+=3,optr+=3){
		*optr+=forward_table[*iptr]/resample;
		/* Red */

		optr[1]+=forward_table[iptr[1]]/resample;
		/* Green */

		optr[2]+=forward_table[iptr[2]]/resample;
		/* Blue */
	}
}

void * mem_alloc(size_t size)
{
	void *ret;

	ret=malloc(size);
	if (size&&!ret){
		fprintf(stderr,"%s: out of memory\n", progname);
		exit(4);
	}
	return ret;
}

/* Needs the line_buffer. return 0=OK, 1=error/end of file.
 * In case of success, allocates d1 and fills in with the
 * loaded data. The d1 contents are photometrically linear RGB.
 * resample is temporal resample. phase is phase of the resample
 * (beginning with 0 */
int load_input(unsigned resample, unsigned phase)
{
	size_t rd;
	unsigned y;
	unsigned short *ptr;

	if (!phase) d1=mem_alloc(3*input_w*input_h*sizeof(*d1));

	if (phase==1) /* Divide by resample */
		for (ptr=d1; ptr<d1+input_w*input_h*3; ptr++){
			*ptr/=resample;
	/* We do not round here because it could case arithmetic overflow. */
		}

	for (y=0;y<input_h;y++){
		rd=fread(line_buffer, 3, input_w, stdin);
		if (rd!=input_w){
			free(d1);
			return 1; /* End of file or error. */
		}
		if (!phase) write_line(y);
		else write_line_div_add(y, resample);
	}
	return 0;
}

void * mem_calloc(size_t size)
{
	void *ret;

	ret=calloc(size, 1);
	if (size&&!ret){
		fprintf(stderr,"%s: out of memory\n", progname);
		exit(4);
	}
	return ret;
}

#define mem_free free

#define CLIP_LOOKUP(x) {if ((x)<0) (x)=0; else if ((x)>255) (x)=255; \
	(x)=forward_table[(x)];}
/* Input: 	input[0]=Y'[0]
 * 		input[1]=Y'[1]
 * 		input[2]=Y'[2]
 * 		input[3]=Y'[3]
 * 		input[4]=Cb
 * 		input[5]=Cr 
 * Attention - changes the input!
 * rc, gc, bc - R'G'B' corresponding to average specified RGB on the input.
 * errors: 4 elements: Yc error, R' error, G' error, B' error
 */
unsigned long perceived_error (unsigned long long *errors, unsigned char *yc
	,unsigned char rc, unsigned char gc, unsigned char bc /* rc,gc,bc are
								 0-255 ! */
	,int *input)
{
	long rx_r, rx_g, rx_b; /* 0-65535 plus must hold some negative
				  values. */
	long rx_r_sum=0, rx_g_sum=0, rx_b_sum=0; /* 0-4*65535 */
	unsigned short y;
	int i;
	unsigned long error_sum=0;
	unsigned long error;

	input[4]-=128; /* Cb-=128 */
	input[5]-=128; /* Cr-=128 */

	for (i=0;i<4;i++){
		input[i]-=16;
		/* Calculate received linear R,G,B components */
		rx_r=SRAL(LUMA_16*input[i]+R_CR_16*input[5]+0x8000,16);
		CLIP_LOOKUP(rx_r);
		rx_g=SRAL(LUMA_16*input[i]+G_CB_16*input[4]+G_CR_16*input[5]+0x8000,16);
		CLIP_LOOKUP(rx_g);
		rx_b=SRAL(LUMA_16*input[i]+B_CB_16*input[4]+0x8000,16);
		CLIP_LOOKUP(rx_b);

		/* Sum into received linear overall R,G,B component
		 * of the 2x2 pixel square */
		rx_r_sum+=rx_r;
		rx_g_sum+=rx_g;
		rx_b_sum+=rx_b;

		y=(LUMI_R_16*rx_r+LUMI_G_16*rx_g+LUMI_B_16*rx_b+0x8000)>>16;
		/* Now y contains luminance, 0-65535 */
		y=backward_table[y];
		/* Now y contains gamma-companded luminance 0-255 */
		error=(int)y-yc[i];
		error=(short)error*error;
		/* Now y contains gamma-companded luminance error */
		if (errors) errors[0]+=error; /* Adds 65025 max. */
		else error_sum+=error;
	}

	/* Normalize the resulting RGB sum to 0-65535 and convert to R'G'B' */
	rx_r_sum=backward_table[(rx_r_sum+2)>>2];
	rx_g_sum=backward_table[(rx_g_sum+2)>>2];
	rx_b_sum=backward_table[(rx_b_sum+2)>>2];

	/* rx_r, rx_g, rx_b now contain gamma-companded cone signals */
	rx_r_sum-=rc;
	rx_g_sum-=gc;
	rx_b_sum-=bc;
	/* rx_r, rx_g, rx_b now contain cone errors */
	if (errors){
		errors[1]+=(short)rx_r_sum*rx_r_sum; /* Adds 65025 max. */
		errors[2]+=(short)rx_g_sum*rx_g_sum; /* Adds 65025 max. */
		errors[3]+=(short)rx_b_sum*rx_b_sum; /* Adds 65025 max. */
	}else
		error_sum+=(short)rx_r_sum*rx_r_sum
			+(short)rx_g_sum*rx_g_sum
			+(short)rx_b_sum*rx_b_sum;

	input[0]+=16;
	input[1]+=16;
	input[2]+=16;
	input[3]+=16;
	input[4]+=128;
	input[5]+=128;
	total_iterations++;
	return error_sum;
}

#define SAVE_LUMA(i) \
	{luma[(i)]=((LUMA_R_16*rc+LUMA_G_16*gc+LUMA_B_16*bc+0x108000)>>16);}
#define SAVE_CB \
	{cb_sum+=CB_R_16*rc+CB_G_16*gc+CB_B_16*bc;}
#define SAVE_CR \
	{cr_sum+=CR_R_16*rc+CR_G_16*gc+CR_B_16*bc;}
#define SAVE_YCBCR(i) {SAVE_LUMA(i) SAVE_CB SAVE_CR}
/* ih - pointer to upper left pixel, contains RGB RGB in unsigned short 0-65535
 * in - pointer to the upper left input pixel's R. Lower left pixel is at
 *      in+3*w.
 * w  - width of image (in pixels) 
 */
void ycbcr_square(
		unsigned char *luma, unsigned char *cb, unsigned char
		*cr,
		unsigned short *in, unsigned w)
{
	unsigned char rc, gc, bc; /* 0-255 */
	unsigned short y[4]; /* Luminance */
	unsigned char yc[4]; /* Gamma-companded luminance */
	unsigned w3=3*w;
	int pei[6]; /* At least 9 bits signed */
	int cb_sum=0, cr_sum=0; /* At least 10 bits signed */

	/* Calculate luminances for all 4 pixels */
	y[0]=(LUMI_R_16*in[0]+LUMI_G_16*in[1]+LUMI_B_16*in[2]+0x8000)>>16;
	y[1]=(LUMI_R_16*in[3]+LUMI_G_16*in[4]+LUMI_B_16*in[5]+0x8000)>>16;
	y[2]=(LUMI_R_16*in[w3]+LUMI_G_16*in[w3+1]+LUMI_B_16*in[w3+2]+0x8000)
		>>16;
	y[3]=(LUMI_R_16*in[w3+3]+LUMI_G_16*in[w3+4]+LUMI_B_16*in[w3+5]+0x8000)
		>>16;

	/* Calculate gamma-companded luminances */
	yc[0]=backward_table[y[0]];
	yc[1]=backward_table[y[1]];
	yc[2]=backward_table[y[2]];
	yc[3]=backward_table[y[3]];

	rc=backward_table[in[0]];
	gc=backward_table[in[1]];
	bc=backward_table[in[2]];
	SAVE_YCBCR(0);
	rc=backward_table[in[3]];
	gc=backward_table[in[4]];
	bc=backward_table[in[5]];
	SAVE_YCBCR(1);
	rc=backward_table[in[w3]];
	gc=backward_table[in[w3+1]];
	bc=backward_table[in[w3+2]];
	SAVE_YCBCR(w);
	rc=backward_table[in[w3+3]];
	gc=backward_table[in[w3+4]];
	bc=backward_table[in[w3+5]];
	SAVE_YCBCR(w+1);

	pei[0]=luma[0];
	pei[1]=luma[1];
	pei[2]=luma[w];
	pei[3]=luma[w+1];

	/* Calculate average original R',G',B' for all 4 pixels together */
	rc=backward_table[((unsigned long)in[0]+in[3]+in[w3]+in[w3+3]+2)>>2];
	gc=backward_table[((unsigned long)in[1]+in[4]+in[w3+1]+in[w3+4]+2)>>2];
	bc=backward_table[((unsigned long)in[2]+in[5]+in[w3+2]+in[w3+5]+2)>>2];

	/* Output the chroma */
	assert(cb_sum>=-0x2020000);
	assert(cr_sum>=-0x2020000);

	*cb=pei[4]=((cb_sum+0x2020000)>>18);
	*cr=pei[5]=((cr_sum+0x2020000)>>18);

	last_error=perceived_error(total_errors
		,yc , rc, gc, bc
		, pei);
	total_quadpixels+=1;
}

/* Only when yc,cb,cr are limited to 0-255, 16-240, 16-240 ! */
/* NOTABLE means the Hyperluma table is not yet calculated and linked. It's
 * replaced with zero just to make the program compile and link. */
#define HYPERLUMA_LOOKUP(yc,cb,cr) \
	hyperluma_table[(yc)+((((cr)-16)*225+((cb)-16))<<8)]

/* in - pointer to upper left pixel, contains RGB RGB in unsigned short 0-65535
 * and in[w*3] constains another RGB RGB 
 * Non-reentrant - stores last result in a static variable and reuses it!
 * Increments total_quadpixels.
 * Adds errors either into total_errors[0-3] (normal operation) or
 * total_errors[4] (reuse). */
void luminaplex_square(
		unsigned char *luma, unsigned char *cb, unsigned char
		*cr,
		unsigned short *in, unsigned w, int algo, int reuse)
{
	unsigned char rc, gc, bc; /* 0-255. Original RGB summed then converted into
				     R'G'B' */
	unsigned rcsum, gcsum, bcsum; /* R'G'B' calculated for each pixel and
					 resulting R'G'B' values summed
					 together. 0-4*255. */
	unsigned char yc[4]; /* Gamma-companded luminance */
	unsigned w3=w*3;
	static int pei[6]; /* Perceived Error Input. As a sideeffect, this
			      contains the results from the last run:
			      pei[0]=luma[0]
			      pei[1]=luma[1]
			      pei[2]=luma[w]
			      pei[3]=luma[w+1]
			      pei[4]=*cb
			      pei[5]=*cr */
	unsigned char *hltable; /* Hyperluma table base */
	unsigned long best_error;
	unsigned long new_error;
	int i;
	unsigned char improved;
	static unsigned short last[12];
	static unsigned char last_vals[6];
	/* Contains: luma[0], luma[1], luma[w], luma[w+1], *cb, *cr */
	static int last_valid; /* 1 means that pei holds the last values,
		last_error holds the last error value and last[12] holds the
		last RGBRGBRGBRGB data and all these data are consistent against
		each other. */

	if (reuse){
		if ((*last==*in) /* This is fast to evaluate and mostly fails */
			&&last_valid /* This passes always except for the
					    first time in case we don't
					  use ALGO_ORDINARY, otherwise it
					  almost always fails. Placing it
					  here to not slow down ALGO_ORDINARY
					  with the memcmp's.
					  */
			&&!memcmp(last+1,in+1,5*sizeof *last)
			&&!memcmp(last+6, in+w3, 6*sizeof *last)
			){

#ifdef DEBUG_HIGHLIGHT_REUSE
			luma[0]=128;
			luma[1]=128;
			luma[w]=128;
			luma[w+1]=128;
			*cb=128;
			*cr=128;
#else
			luma[0]=last_vals[0];
			luma[1]=last_vals[1];
			luma[w]=last_vals[2];
			luma[w+1]=last_vals[3];
			*cb=last_vals[4];
			*cr=last_vals[5];
#endif
			total_errors[4]+=last_error;
			total_quadpixels++;
			/* Repetition */
			repetitions++;
			return;
		}else if (reuse){
			memcpy(last, in, 6*sizeof *last);
			memcpy(last+6, in+w3, 6*sizeof *last);
			last_valid=0; /* For a while it's gonna be invalid.
					     If we call ALGO_ORDINARY, then it
					     stays invalid. Therefore it has
					     to be invalidated. */
			non_repetitions++; /* This one apparently won't be
					      copied from the previous one ;-)
					    */
		}
	}
	if (algo==ALGO_ORDINARY){
		ycbcr_square(luma, cb, cr, in, w);
		goto store_for_future;
	}
	/* Calculate gamma-companded luminances for all 4 pixels */
	yc[0]=backward_table[(LUMI_R_16*in[0]+LUMI_G_16*in[1]+LUMI_B_16*in[2]+0x8000)>>16];
	yc[1]=backward_table[(LUMI_R_16*in[3]+LUMI_G_16*in[4]+LUMI_B_16*in[5]+0x8000)>>16];
	yc[2]=backward_table[(LUMI_R_16*in[w3]+LUMI_G_16*in[w3+1]+LUMI_B_16*in[w3+2]+0x8000)
		>>16];
	yc[3]=backward_table[(LUMI_R_16*in[w3+3]+LUMI_G_16*in[w3+4]+LUMI_B_16*in[w3+5]+0x8000)
		>>16];

	/* Average the input RGB for all 4 pixels and then convert it into
	 * R'G'B'. */
	rc=backward_table[((unsigned long)in[0]+in[3]+in[w3]+in[w3+3]+2)>>2];
	gc=backward_table[((unsigned long)in[1]+in[4]+in[w3+1]+in[w3+4]+2)>>2];
	bc=backward_table[((unsigned long)in[2]+in[5]+in[w3+2]+in[w3+5]+2)>>2];

	/* Convert the input RGB into R'G'B' for each pixel and then sum the
	 * results. Max. result: 4*255 */
	rcsum=backward_table[in[0]]+backward_table[in[3]]+backward_table[in[w3]]+backward_table[in[w3+3]];
	gcsum=backward_table[in[1]]+backward_table[in[4]]+backward_table[in[w3+1]]+backward_table[in[w3+4]];
	bcsum=backward_table[in[2]]+backward_table[in[5]]+backward_table[in[w3+2]]+backward_table[in[w3+5]];

	/* rc, gc, bc now range 0-65535. */
	/* rcsum, gcsum, bcsum now range 0-255*4 */

	/* Calculate the chroma */
	if (algo==ALGO_HYPERLUMA_1){
	pei[4]=((CB_R_16*rcsum+CB_G_16*gcsum+CB_B_16*bcsum+0x2020000)>>18);
	pei[5]=((CR_R_16*rcsum+CR_G_16*gcsum+CR_B_16*bcsum+0x2020000)>>18);
	}else{
		/* ALGO_HYPERLUMA_2 or ALGO_LUMINAPLEX */
	pei[4]=((CB_R_16*rc+CB_G_16*gc+CB_B_16*bc+0x808000)>>16);
	pei[5]=((CR_R_16*rc+CR_G_16*gc+CR_B_16*bc+0x808000)>>16);
	}

	/* Calculate the 4 lumas */
	hltable=hyperluma_table+(((pei[5]-16)*225UL+(pei[4]-16))<<8);
	pei[0]=hltable[yc[0]];
	pei[1]=hltable[yc[1]];
	pei[2]=hltable[yc[2]];
	pei[3]=hltable[yc[3]];

	if (algo==ALGO_LUMINAPLEX){
		best_error=perceived_error(NULL, yc, rc, gc, bc, pei);
		do
		{
			unsigned char loc_imp; /* Local improved */

			improved=0;
			for (i=0;i<4;i++){
				loc_imp=0;
				do {
					pei[i]++;
					if (pei[i]>235) break;
					new_error=perceived_error(NULL, yc, rc,
							gc, bc, pei);
					if (new_error>best_error) break;
					if (new_error<best_error){
						improved=1;
						loc_imp=1;
					}
					best_error=new_error;
				}while(1);
				pei[i]--;

				if (!loc_imp){
					do {
						pei[i]--;
						if (pei[i]<16) break;
						new_error=perceived_error(NULL, yc, rc,
								gc, bc, pei);
						if (new_error>best_error) break;
						if (new_error<best_error) improved=1;
						best_error=new_error;
					}while(1);
					pei[i]++;
				}

			}
			for (i=4;i<6;i++){
				loc_imp=0;
				do {
					pei[i]++;
					if (pei[i]>240) break;
					new_error=perceived_error(NULL, yc, rc,
							gc, bc, pei);
					if (new_error>best_error) break;
					if (new_error<best_error){
						improved=1;
						loc_imp=1;
					}
					best_error=new_error;
				}while(1);
				pei[i]--;

				if (!loc_imp){
					do {
						pei[i]--;
						if (pei[i]<16) break;
						new_error=perceived_error(NULL, yc, rc,
								gc, bc, pei);
						if (new_error>best_error) break;
						if (new_error<best_error) improved=1;
						best_error=new_error;
					}while(1);
					pei[i]++;
				}

			}
		}while(improved);
	}
	luma[0]=pei[0];
	luma[1]=pei[1];
	luma[w]=pei[2];
	luma[w+1]=pei[3];
	*cb=pei[4];
	*cr=pei[5];
	/* Prepare input for the perceived error routine */
	last_error=perceived_error(total_errors, yc, rc, gc, bc, pei);
	total_quadpixels+=1;
store_for_future:
	last_vals[0]=luma[0];
	last_vals[1]=luma[1];
	last_vals[2]=luma[w];
	last_vals[3]=luma[w+1];
	last_vals[4]=*cb;
	last_vals[5]=*cr;
	last_valid=1;
	return;
}

void write_test_line(unsigned short *ptr, unsigned short r,
		unsigned short g, unsigned short b, unsigned pixels)
{
	unsigned short *endptr;

	for (endptr=ptr+3*pixels; ptr<endptr; ptr+=3){
		*ptr=r;
		ptr[1]=g;
		ptr[2]=b;
	}
}

void export_rgb_16(unsigned short *inbuf, unsigned w, unsigned h)
{
	unsigned char *outbuf=mem_alloc(3*w);
	unsigned char *outbufptr;
	unsigned short *line_begin, *inbufptr;

	line_begin=inbuf+3*w*(h-1);
	for (;h;h--){
		outbufptr=outbuf;
		inbufptr=line_begin;
		while(outbufptr<outbuf+3*w)
			*outbufptr++=backward_table[*inbufptr++];
		my_fwrite(outbuf, 3, w, stdout);
		line_begin-=3*w;
	}

	mem_free(outbuf);
}

#define CSTEP(x) {*table++=color; color=(color+x)&7;}

void make_colors(unsigned char *table)
{
	int i;
	unsigned color=0;

	for (i=8; i; i--){
		CSTEP(1);
		CSTEP(2);
		CSTEP(4);
	}
	for (i=8; i; i--){
		CSTEP(3);
	}
	*table=color;

}

void sim_quadpixel(unsigned short *top, unsigned short *bottom, int algo)
{
	unsigned short in[12];
	unsigned char out[6];
	int i;
	long rx_r, rx_g, rx_b;
	int y, cb, cr;

	if (algo>=4) return;
	memcpy(in, top, 6*sizeof(*in));
	memcpy(in+6, bottom, 6*sizeof(*in));
	luminaplex_square(out, out+4, out+5, in, 2, algo,1);

	cb=out[4]-128;
	cr=out[5]-128;

	for (i=0;i<4; i++){
		
		y=out[i]-16;
		/* Calculate received linear R,G,B components */
		rx_r=SRAL(LUMA_16*y+R_CR_16*cr+0x8000,16);
		CLIP_LOOKUP(rx_r);
		rx_g=SRAL(LUMA_16*y+G_CB_16*cb+G_CR_16*cr+0x8000,16);
		CLIP_LOOKUP(rx_g);
		rx_b=SRAL(LUMA_16*y+B_CB_16*cb+0x8000,16);
		CLIP_LOOKUP(rx_b);
		in[i*3]=rx_r;
		in[i*3+1]=rx_g;
		in[i*3+2]=rx_b;
	}
	memcpy(top, in, 6*sizeof(*in));
	memcpy(bottom, in+6, 6*sizeof(*in));
}

#define TEST_WIDTH 320 /* quadpixels */
#define TEST_HEIGHT 240 /* quadpixels */
void testcard(void)
{
	unsigned w=TEST_WIDTH<<1;
	unsigned h=TEST_HEIGHT<<1;
	unsigned short *rgb=mem_alloc(w*h*3*sizeof(*rgb));
	unsigned short *ptr;
	unsigned bar_index; 
	unsigned y,x;
	unsigned char colors[33];
	unsigned white,black;
	
	white=desaturate?forward_table[0xd0]:0xffff;
	black=desaturate?forward_table[0x30]:0;

	make_colors(colors);

	for (ptr=rgb,y=0;y<h; ptr+=w*3,y++){
		bar_index=colors[(((y+1)&-2)*33)/(h|1)];
		write_test_line(ptr
				,(bar_index&2)?white:black
				,(bar_index&4)?white:black
				,(bar_index&1)?white:black
				,w);
	}
	for (y=0;y+1<h;y+=2){
		for (x=0,ptr=rgb+3*w*y;x+1<w;ptr+=6,x+=2)
			sim_quadpixel(ptr, ptr+3*w, (x*5)/w);

	}
	export_rgb_16(rgb, w,h);
}

void luminaplex_line(unsigned short *src, unsigned char *luma,
		unsigned char *cb, unsigned char *cr, unsigned w)
{
	unsigned char *cbend;
	unsigned short in_tmp[12];
	unsigned char luma_tmp[4];


	for (cbend=cb+(w>>1);cb<cbend;cb++){
		luminaplex_square(luma, cb, cr, src, w, ALGO, 1);
		src+=6;
		luma+=2;
		cr++;
	}
	if (w&1){
		memcpy(in_tmp, src, 3*sizeof(*in_tmp));
		memcpy(in_tmp+3, src, 3*sizeof(*in_tmp));
		memcpy(in_tmp+6, src+3*w, 3*sizeof(*in_tmp));
		memcpy(in_tmp+9, src+3*w, 3*sizeof(*in_tmp));
		luminaplex_square(luma_tmp, cb, cr, in_tmp, 2, ALGO, 1);
		*luma=luma_tmp[0];
		luma[w]=luma_tmp[2];
	}
}

void luminaplex_lastline(unsigned short *src, unsigned char *luma,
		unsigned char *cb, unsigned char *cr, unsigned w)
{
	unsigned char *cbend;
	unsigned short in_tmp[12];
	unsigned char luma_tmp[4];


	for (cbend=cb+(w>>1);cb<cbend;cb++){
		memcpy(in_tmp, src, 6*sizeof(*in_tmp));
		memcpy(in_tmp+6, src, 6*sizeof(*in_tmp));
		luminaplex_square(luma_tmp, cb, cr, in_tmp, 2, ALGO, 1);
		luma[0]=luma_tmp[0];
		luma[1]=luma_tmp[1];
		src+=6;
		luma+=2;
		cr++;
	}
	if (w&1){
		memcpy(in_tmp, src, 3*sizeof(*in_tmp));
		memcpy(in_tmp+3, src, 3*sizeof(*in_tmp));
		memcpy(in_tmp+6, src, 3*sizeof(*in_tmp));
		memcpy(in_tmp+9, src, 3*sizeof(*in_tmp));
		luminaplex_square(luma_tmp, cb, cr, in_tmp, 2, ALGO, 1);
		*luma=luma_tmp[0];
	}
}

void luminaplex_array(unsigned char *dest, unsigned short *src, unsigned w,
		unsigned h)
{
	unsigned char *luma, *chroma;
	unsigned chromaoff; /* Offset between cb and cr */
	int i=0;

	luma=dest;
	chroma=luma+w*h;
	chromaoff=((w+1)>>1)*((h+1)>>1);
	for (;luma<dest+w*(h&-2);luma+=(w<<1)){
		luminaplex_line(src, luma, chroma, chroma+chromaoff,
				w);
		i+=2;
		src+=6*w;	
		chroma+=(w+1)>>1;
	}
	if (h&1) luminaplex_lastline(src, luma, chroma
			,chroma+chromaoff, w);
}

void dump_luminaplex_frame(unsigned short *in, unsigned w, unsigned h)
{
	unsigned char *out; /* Luma-chroma-chroma planes */
	unsigned outlen; /* Length of the luma and the 2 chroma planes together
			  */

	outlen=w*h+(((w+1)>>1)*((h+1)&-2));
	out=mem_alloc(outlen);
	luminaplex_array(out, in, w, h);
	write_frame_header();
	my_fwrite(out, outlen, 1, stdout);
	mem_free(in);
	mem_free(out);
}

void make_gamma_tables(void){
	int a;
	float f;

	for (a=0;a<0x100;a++){
		forward_table[a]=floor(0xffff*pow(a/255.0, 1/input_gamma)+0.5);
	}
	for (a=0;a<0x10000;a++){
		f=pow(a/65535.0, 1/TV_GAMMA);
		backward_table_16[a]=floor(0xffff*f+0.5);
		backward_table[a]=floor(0xff*f+0.5);
	}
}

#define overalloc()							\
do {									\
	fprintf(stderr,\
	"%s: ERROR: attempting to allocate too large block at %s:%d", \
	progname, __FILE__, __LINE__);\
	exit(5);						\
} while (1)	/* while (1) is not a typo --- it's here to allow a
	compiler that doesn't know that exit doesn't return to do better
	optimizations */

/* We assume unsigned holds at least 32 bits */
inline static void bias_buf_color(unsigned *col_buf, int n, unsigned half)
{
	for (;n;n--){
		*col_buf=half;
		col_buf[1]=half;
		col_buf[2]=half;
		col_buf+=3;
	}
}
		
/* line_skip is in pixels. The column contains the whole pixels (R G B)
 * We assume unsigned short holds at least 16 bits. */
inline static void add_col_color(unsigned *col_buf, unsigned short *ptr
	, int line_skip, int n, unsigned weight)
{
	for (;n;n--){
		*col_buf+=weight*(*ptr);
		col_buf[1]+=weight*ptr[1];
		col_buf[2]+=weight*ptr[2];
		ptr+=line_skip;
		col_buf+=3;
	}
}

/* line skip is in pixels. Pixel is 3*unsigned short */
/* We assume unsigned holds at least 32 bits */
/* We assume unsigned short holds at least 16 bits. */
inline static void emit_and_bias_col_color(unsigned *col_buf
	, unsigned short *out, int line_skip, int n, unsigned weight)
{
	unsigned half=weight>>1;

	for (;n;n--){
		*out=(*col_buf)/weight;
		*col_buf=half;
		out[1]=col_buf[1]/weight;
		col_buf[1]=half;
		/* The following line is an enemy of the State and will be
		 * prosecuted according to the Constitution of The United States
		 * Cap. 20/3 ix. Sel. Bill 12/1920
		 * Moses 12/20 Erizea farizea 2:2:1:14
		 */
		out[2]=col_buf[2]/weight;
		col_buf[2]=half;
		out+=line_skip;
		col_buf+=3;
	}
}
		
/* For enlargement only -- does linear filtering
 * Frees input and allocates output.
 * We assume unsigned holds at least 32 bits
 */
static inline void enlarge_color_horizontal(unsigned short *ina, int ix, int y,
	unsigned short ** outa, int ox)
{
	unsigned *col_buf;
	int total,a,out_pos,in_pos,in_begin,in_end;
	unsigned half=(ox-1)>>1;
	unsigned skip=3*ix;
	unsigned oskip=3*ox;
	unsigned short *out, *in;

	if (ix==ox){
		*outa=ina;
		return;
	}
	if (ox && (unsigned)ox * (unsigned)y / (unsigned)ox != (unsigned)y) overalloc();
	if ((unsigned)ox * (unsigned)y > UINT_MAX / 3 / sizeof(*out)) overalloc();
	out=mem_alloc(sizeof(*out)*3*ox*y);
	*outa=out;
	in=ina;
	if (ix==1){
		for (;y;y--,in+=3) for (a=ox;a;a--,out+=3){
			*out=*in;
			out[1]=in[1];
			out[2]=in[2];
		}
		mem_free(ina);
		return;
	}
	total=(ix-1)*(ox-1);
	if ((unsigned)y > UINT_MAX / 3 / sizeof(*col_buf)) overalloc();
	col_buf=mem_alloc(y*3*sizeof(*col_buf));
	bias_buf_color(col_buf,y,half);
	out_pos=0;
	in_pos=0;
	again:
	in_begin=in_pos;
	in_end=in_pos+ox-1;
	/* Banzai Pipeline */
	add_col_color(col_buf,in,skip,y
		,in_end-out_pos);
	add_col_color(col_buf,in+3,skip,y
		,out_pos-in_begin);
	emit_and_bias_col_color(col_buf,out,oskip,y,ox-1);
	out+=3;
	out_pos+=ix-1;
	if (out_pos>in_end){
		in_pos=in_end;
		in+=3;
	}
	if (out_pos>total){
		mem_free(col_buf);
		mem_free(ina);
		return;
	}
	goto again;
}

/* n is in pixels. pixel is 3 unsigned shorts in series */
 /* We assume unsigned short holds at least 16 bits. */
inline static void add_row_color(unsigned *row_buf, unsigned short *ptr, int n, unsigned weight)
{
	for (;n;n--){
		*row_buf+=weight**ptr;
		row_buf[1]+=weight*ptr[1];
		row_buf[2]+=weight*ptr[2];
		ptr+=3;
		row_buf+=3;
	}
}

/* n is in pixels. pixel is 3 unsigned shorts in series. */
/* We assume unsigned holds at least 32 bits */
/* We assume unsigned short holds at least 16 bits. */
inline static void emit_and_bias_row_color(unsigned *row_buf, unsigned short
		*out, int n, unsigned weight)
{
	unsigned half=weight>>1;

	for (;n;n--){
		*out=*row_buf/weight;
		*row_buf=half;
		out[1]=row_buf[1]/weight;
		row_buf[1]=half;
		out[2]=row_buf[2]/weight;
		row_buf[2]=half;
		out+=3;
		row_buf+=3;
	}
}
		
/* For magnification only. Does linear filtering */
/* We assume unsigned holds at least 32 bits */
inline static void enlarge_color_vertical(unsigned short *ina, int x, int iy,
	unsigned short **outa ,int oy)
{
	unsigned *row_buf;
	int total,out_pos,in_pos,in_begin,in_end;
	int half=(oy-1)>>1;
	unsigned short *out, *in;

	if (iy==oy){
		*outa=ina;
		return;
	}
	/* Pacific Ocean Park */
	if (x && (unsigned)x * (unsigned)oy / (unsigned)x != (unsigned)oy) overalloc();
	if ((unsigned)x * (unsigned)oy > UINT_MAX / 3 / sizeof(*out)) overalloc();
	out=mem_alloc(sizeof(*out)*3*oy*x);
	*outa=out;
	in=ina;
	if (iy==1){
		for (;oy;oy--){
	       		memcpy(out,in,3*x*sizeof(*out));
	       		out+=3*x;
		}
		mem_free(ina);
		return;
	}
	total=(iy-1)*(oy-1);
	if ((unsigned)x > UINT_MAX / 3 / sizeof(*row_buf)) overalloc();
	row_buf=mem_alloc(x*3*sizeof(*row_buf));
	bias_buf_color(row_buf,x,half);
	out_pos=0;
	in_pos=0;
	again:
	in_begin=in_pos;
	in_end=in_pos+oy-1;
	add_row_color(row_buf,in,x
		,in_end-out_pos);
	add_row_color(row_buf,in+3*x,x
		,out_pos-in_begin);
	emit_and_bias_row_color(row_buf,out,x,oy-1);
	out+=3*x;
	out_pos+=iy-1;
	if (out_pos>in_end){
		in_pos=in_end;
		in+=3*x;
	}
	if (out_pos>total){
		mem_free(ina);
		mem_free(row_buf);
		return;
	}
	goto again;
	
}	

/* Works for both enlarging and diminishing. Linear resample, no low pass.
 * Does only one color component.
 * Frees ina and allocates outa.
 * If ox*3<=ix, and display_optimize, performs optimization for LCD.
 */
inline static void scale_color_horizontal(unsigned short *ina, int ix, int y,
		unsigned short **outa, int ox)
{
	unsigned *col_buf;
	int total=ix*ox;
	int out_pos,in_pos,in_begin,in_end,out_end;
	unsigned skip=3*ix;
	unsigned oskip=3*ox;
	unsigned short *in, *out;

	if (ix==ox){
		*outa=ina;
		return;
	}
	if (ix<ox){
		enlarge_color_horizontal(ina,ix,y,outa,ox);
		return;
	}else if (ix==ox){
		*outa=ina;
		return;
	}
	if (ox && (unsigned)ox * (unsigned)y / (unsigned)ox != (unsigned)y) overalloc();
	if ((unsigned)ox * (unsigned)y > UINT_MAX / 3 / sizeof(*out)) overalloc();
	out=mem_alloc(sizeof(*out)*3*ox*y);
	*outa=out;
	in=ina;
	if ((unsigned)y > UINT_MAX / 3 / sizeof(*col_buf)) overalloc();
	col_buf=mem_alloc(y*3*sizeof(*col_buf));
	bias_buf_color(col_buf,y,ix>>1);
	out_pos=0;
	in_pos=0;
	again:
	in_begin=in_pos;
	in_end=in_pos+ox;
	out_end=out_pos+ix;
	if (in_begin<out_pos)in_begin=out_pos;
	if (in_end>out_end)in_end=out_end;
	add_col_color(col_buf,in,skip,y,in_end-in_begin);
	in_end=in_pos+ox;
	if (out_end>=in_end){
		in_pos=in_end;
		in+=3;
	}
	if (out_end<=in_end){
			emit_and_bias_col_color(col_buf,out,oskip,y,ix);
			out_pos=out_pos+ix;
			out+=3;
	}
	if (out_pos==total) {
		mem_free(ina);
		mem_free(col_buf);
		return;
	}
	goto again;
}

/* Both enlarges and diminishes. Linear filtering. Sizes are
   in pixels. Sizes are not in bytes. 1 pixel=3 unsigned shorts.
   We assume unsigned short can hold at least 16 bits.
   We assume unsigned holds at least 32 bits.
 */
inline static void scale_color_vertical(unsigned short *ina, int x, int iy
	,unsigned short **outa, int oy)
{
	unsigned *row_buf;
	int total=iy*oy;
	int out_pos,in_pos,in_begin,in_end,out_end;
	unsigned short *in, *out;

	if (iy==oy){
		*outa=ina;
		return;
	}
	if (iy<oy){
		enlarge_color_vertical(ina,x,iy,outa,oy);
		return;
	}
	if (x && (unsigned)x * (unsigned)oy / (unsigned)x != (unsigned)oy) overalloc();
	if ((unsigned)x * (unsigned)oy > UINT_MAX / 3 / sizeof(*out)) overalloc();
	out=mem_alloc(sizeof(*out)*3*oy*x);
	*outa=out;
	in=ina;
	if ((unsigned)x > UINT_MAX / 3 / sizeof(*row_buf)) overalloc();
	row_buf=mem_alloc(x*3*sizeof(*row_buf));
	bias_buf_color(row_buf,x,iy>>1);
	out_pos=0;
	in_pos=0;
	again:
	in_begin=in_pos;
	in_end=in_pos+oy;
	out_end=out_pos+iy;
	if (in_begin<out_pos)in_begin=out_pos;
	if (in_end>out_end)in_end=out_end;
	add_row_color(row_buf,in,x,in_end-in_begin);
	in_end=in_pos+oy;
	if (out_end>=in_end){
		in_pos=in_end;
		in+=3*x;
	}
	if (out_end<=in_end){
			emit_and_bias_row_color(row_buf,out,x,iy);
			out_pos=out_pos+iy;
			out+=3*x;
	}
	if (out_pos==total){
		mem_free(ina);
		mem_free(row_buf);
		return;
	}
	goto again;
}


/* Scales color 48-bits-per-pixel bitmap. Both enlarges and diminishes. Uses
 * either low pass or bilinear filtering. The memory organization for both
 * input and output are red, green, blue. All three of them are unsigned shorts 0-65535.
 * Allocates output and frees input
 * We assume unsigned short holds at least 16 bits.
 */
void scale_color(unsigned short *in, int ix, int iy, unsigned short **out,
	int ox, int oy)
{
	unsigned short *intermediate_buffer;

	if (!ix||!iy){
		if (in) mem_free(in);
		if (ox && (unsigned)ox * (unsigned)oy / (unsigned)ox != (unsigned)oy) overalloc();
		if ((unsigned)ox * (unsigned)oy > UINT_MAX / 3 / sizeof(**out)) overalloc();
		*out=mem_calloc(ox*oy*sizeof(**out)*3);
		return;
	}
	if (ix*oy<ox*iy){
		scale_color_vertical(in,ix,iy,&intermediate_buffer,oy);
		scale_color_horizontal(intermediate_buffer,ix,oy,out,ox);
	}else{
		scale_color_horizontal(in,ix,iy,&intermediate_buffer,ox);
		scale_color_vertical(intermediate_buffer,ox,iy,out,oy);
	}
}

/* Converts photometrically linear 48bpp RGB to gamma=1/TV_GAMMA 24bpp plane-separated
 * R'G'B'. Allocates d3 and frees d2. */
void d2_to_d3(void)
{
	unsigned offset=output_w*output_h; /* Offset between green and red write pointer and
				also between blue and green write pointer */
	unsigned char *wptr; /* Write pointer */
	unsigned short *rptr; /* Read pointer */

	d3=mem_alloc(offset*3);
	if (ppm){
		for (rptr=d2,wptr=d3; wptr<d3+offset*3; wptr+=3, rptr+=3){
			*wptr=backward_table[*rptr]; /* Red */
			wptr[1]=backward_table[rptr[1]]; /* Green */
			wptr[2]=backward_table[rptr[2]]; /* Blue */
		}
	}else{
		for (rptr=d2,wptr=d3; wptr<d3+offset; wptr++, rptr+=3){
			*wptr=backward_table[*rptr]; /* Red */
			wptr[offset]=backward_table[rptr[1]]; /* Green */
			wptr[offset<<1]=backward_table[rptr[2]]; /* Blue */
		}
	}
	mem_free(d2);
}

/* Doesn't touch any allocation. */
void dump_ppm_frame(unsigned char *d)
{
	printf("P6 %u %u 255\n", output_w, output_h);
	my_fwrite(d, output_w*3, output_h, stdout);
}

void store_rle_block(unsigned long *position, unsigned char c, unsigned count)
{
	static int last_value;

	if (*position+count>sizeof hyperluma_table/sizeof *hyperluma_table){
		fprintf(stderr,"%s: Hyperluma file corrupted - too much data\n",
				progname);
		exit(1);
	}
	for (;count;count--){
		last_value=(last_value+c)&255;
		hyperluma_table[(*position)++]=last_value;
	}
}

#ifndef NOTABLE
void uncompress_hyperluma_table(void)
{
	int c;
	unsigned long hyperluma_table_position=0;
	unsigned char *compressed_data_pointer=hyperluma_table_compressed;
	

	/* bit 7: 0: lowest 2 bits value, bits 2-6 repetitions
	 *        1: bits 0-6 repetitions, following byte value
	 */
	while(compressed_data_pointer<hyperluma_table_compressed_end){
		c=*compressed_data_pointer++;
		if (c==EOF){
			fprintf(stderr,"%s: Compressed hyperluma table corrupted "
					"- ran out of "
					"input data\n", progname);
			exit(1);
		}
		if (c&0x80){
			if (compressed_data_pointer
				>=hyperluma_table_compressed_end) break;
			store_rle_block(&hyperluma_table_position
				,*compressed_data_pointer++ , c&~0x80);
		}else{
			store_rle_block(&hyperluma_table_position,c&3, c>>2);
		}
	}
	if (hyperluma_table_position!=sizeof hyperluma_table/sizeof
			*hyperluma_table){
		fprintf(stderr,"%s: Compressed hyperluma data corrupted - too little "
				"data - %lu bytes only\n", progname
				,hyperluma_table_position);
		exit(1);
	}
}
#endif

void square_block(unsigned short *d, unsigned w, unsigned square,
		unsigned x, unsigned y, unsigned xs, unsigned ys,
		unsigned r, unsigned g, unsigned b)
{
	unsigned short *lineptr, *ptr;

	x*=square;
	y*=square;
	xs*=square;
	ys*=square;

	for (lineptr=d+3*(y*w+x);lineptr<d+3*((y+ys)*w+x); lineptr+=3*w)
		for (ptr=lineptr; ptr<lineptr+3*xs; ptr+=3){
			ptr[0]=r;
			ptr[1]=g;
			ptr[2]=b;
		}
	
}

#define WIDTH_TO_SQUARE 70
void hyperluma_logo (unsigned short *d, unsigned w, unsigned h)
{
	int square=(w+WIDTH_TO_SQUARE-1)/WIDTH_TO_SQUARE;
	int xpos=3, ypos=3;

	if ((xpos+8)*square>w||(ypos+6)*square>h) return;
	/* Logo wouldn't fit */

	square_block(d, w, square, xpos-1,ypos-1,7,3,0,0,0);
	square_block(d, w, square, xpos-1,ypos+2,9,4,0,0,0);
	square_block(d, w, square, xpos,ypos,2,5,65535, 0, 65535); /* H */
	square_block(d, w, square, xpos+2,ypos+2,1,2,65535, 0, 65535); /* H */
	square_block(d, w, square, xpos+3,ypos,1,5,65535, 0, 65535); /* H */
	square_block(d, w, square, xpos+4,ypos,1,5,0, 65535, 0); /* L */
	square_block(d, w, square, xpos+5,ypos+3,2,2,0, 65535, 0); /* L */
}

void luminaplex_logo (unsigned short *d, unsigned w, unsigned h)
{
	int square=(w+WIDTH_TO_SQUARE-1)/WIDTH_TO_SQUARE;
	int xpos=3, ypos=3;

	if ((xpos+8)*square>w||(ypos+6)*square>h) return;
	/* Logo wouldn't fit */

	square_block(d, w, square, xpos-1,ypos-1,8,5,0,0,0);
	square_block(d, w, square, xpos-1,ypos+4,6,2,0,0,0);
	square_block(d, w, square, xpos,ypos,1,5,0xffff,0xffff, 0xffff); /* L */
	square_block(d, w, square, xpos+1,ypos+4
			,2,1,0xffff, 0xffff, 0xffff); /* L */
	square_block(d, w, square, xpos+3,ypos,1,5, 65535, 0, 0); /* P */
	square_block(d, w, square, xpos+4,ypos,2,1, 65535, 0, 0); /* P */
	square_block(d, w, square, xpos+5,ypos+1,1,1, 65535, 0, 0); /* P */
	square_block(d, w, square, xpos+4,ypos+2,2,1, 65535, 0, 0); /* P */
}

unsigned long long total_error(void)
{
	return total_errors[0]+total_errors[1]
		+total_errors[2]+total_errors[3]+total_errors[4];
}

/* Returns iterations per pixel */
double iter(void)
{
	return (double)total_iterations/total_quadpixels/4;
}

/* Calculates the RMS error */
double rms_err(double *db)
{
	double rms;

	rms=sqrt((float)total_error()/total_quadpixels/7);
	if (db) *db=20*log10f(255.0/2/rms);
	return rms;
}

/* Return microseconds per pixel */
double usec(void)
{
	struct itimerval current;
	long sec_diff;
	long usec_diff;
	double rv;

	getitimer(ITIMER_USED, &current);
	sec_diff=current.it_interval.tv_sec-current.it_value.tv_sec;
	usec_diff=current.it_interval.tv_usec-current.it_value.tv_usec;
	rv=(sec_diff*(double)1000000
		+(usec_diff))/total_quadpixels/4;
	return rv;
}

void print_errors(void)
{
	double db;
	

	printf("errrors Yc %08llx R'G'B' %07llx %07llx %07llx"
			" total err. %09llx, %05llx quadpixels\n"
			,total_errors[0]
			,total_errors[1]
			,total_errors[2]
			,total_errors[3]
			,total_error()
			,total_quadpixels);
	printf("RMS error=%.3f LSB (how much to add to each 8-bit R',G',B'"
			" channel to get the same simulated error perception)\n"
			,rms_err(&db));
	printf("SNR=%.2f dB, %.2lf usec/pixel, %.2lf iter/pixel\n"
			,db,usec(),iter());
}

void print_result(unsigned char *name, unsigned char *in)
{
	printf("%12s ",name);
	printf("Y' %02x,%02x,%02x,%02x ",in[0], in[1], in[2], in[3]);
	printf("Cb %02x Cr %02x ",in[4], in[5]);
	printf("errrors Yc %04llx R'G'B' %04llx %04llx %04llx total %04llx\n"
			,total_errors[0]
			,total_errors[1]
			,total_errors[2]
			,total_errors[3]
			,total_error());
}

void clear_stats(void)
{
	struct itimerval start_itimer={
		.it_interval={
			.tv_sec=LONG_MAX,
			.tv_usec=0
		},
		.it_value={
			.tv_sec=LONG_MAX,
			.tv_usec=0
		}
	};

	memset(total_errors, 0, sizeof(total_errors));
	total_quadpixels=0;
	total_iterations=0;
	repetitions=0;
	non_repetitions=0;
retry:
	if (setitimer(ITIMER_USED, &start_itimer, NULL)<0){
		if (errno==EINVAL){
			/* Time too large to handle, let's try half */
			fprintf(stderr,"Oh! The interval timer cannot handle "
				"%ld seconds! Trying half...\n"
				, start_itimer.it_value.tv_sec);
			start_itimer.it_value.tv_sec>>=1;
			start_itimer.it_interval.tv_sec>>=1;
			goto retry;
		}
		fprintf(stderr,"%s: failed setting interval timer: "
			, progname);
		perror(NULL);
		exit(1);
	}

}

void test_fn(unsigned short *in)
{
	unsigned char out[6];


	printf("RGB %04x %04x %04x   %04x %04x %04x\n"
			"    %04x %04x %04x   %04x %04x %04x\n"
			,in[0], in[1], in[2], in[3], in[4], in[5]
			,in[6], in[7], in[8], in[9], in[10], in[11]);
	clear_stats();
	ycbcr_square(out, out+4, out+5, in, 2);
	print_result((unsigned char *)"ordinary", out);
	clear_stats();
	luminaplex_square(out, out+4, out+5, in, 2,1,0);
	print_result((unsigned char *)"Hyperluma 1", out);
	clear_stats();
	luminaplex_square(out, out+4, out+5, in, 2,2,0);
	print_result((unsigned char *)"Hyperluma 2", out);
	clear_stats();
	luminaplex_square(out, out+4, out+5, in, 2,3,0);
	print_result((unsigned char *)"Luminaplex", out);
	putchar('\n');
}

int selftest_error(int i)
{
	fprintf(stderr,"%s: selftest failed with %u!\n"
		, progname, i);
	exit(1);
}

void selftest(void)
{
	if (LUMI_R_16+LUMI_G_16+LUMI_B_16!=65536) selftest_error(1);
	if (CB_R_16+CB_G_16+CB_B_16!=0) selftest_error(2);
	if (CR_R_16+CR_G_16+CR_B_16!=0) selftest_error(3);
}

/* 0= passed, 1=failed */
int algo_stats(void)
{
	unsigned short rc, gc, bc;
	unsigned short in[12];
	unsigned short *ptr;
	unsigned char scratchpad[6];
	int i;

	printf("--- COLOR BARS ---\n\n");
	for (i=0;i<8;i++){
		bc=-(i&1);
		rc=-((i>>1)&1);
		gc=-((i>>2)&1);

		in[0]=in[3]=in[6]=in[9]=rc;
		in[1]=in[4]=in[7]=in[10]=gc;
		in[2]=in[5]=in[8]=in[11]=bc;

		test_fn(in);

	}

	printf("--- RANDOM UNIFORM QUADPIXELS ---\n\n");
	for (i=0;i<8;i++){
		rc=random();
		gc=random();
		bc=random();

		in[0]=in[3]=in[6]=in[9]=rc;
		in[1]=in[4]=in[7]=in[10]=gc;
		in[2]=in[5]=in[8]=in[11]=bc;

		test_fn(in);
	}

	printf("--- RANDOM NOISE ---\n\n");
	for (i=0;i<80;i++){
		for (ptr=in; ptr<in+12; ptr++)
			*ptr=random();
		test_fn(in);
	}

	printf("--- POISONOUS QUADPIXEL ---\n\n");
	in[0]=0x0f76;
	in[1]=0x255a;
	in[2]=0xf92e;
	in[3]=0x7263;
	in[4]=0xc233;
	in[5]=0xd79f;
	in[6]=0xc4c9;
	in[7]=0x079a;
	in[8]=0xfb66;
	in[9]=0x5d32;
	in[10]=0x500d;
	in[11]=0xd7b7;
	test_fn(in);

	printf("--- GREEN-MAGENTA EDGE ---\n\n");
	in[0]=in[2]=in[4]=in[6]=in[8]=in[10]=0;
	in[1]=in[3]=in[5]=in[7]=in[9]=in[11]=0xffff;
	test_fn(in);
	
	printf("--- BLACK ---\n\n");
	memset(in, 0, sizeof(in));
	test_fn(in);
	
	printf("--- GREEN-BLACK EDGE ---\n\n");
	memset(in, 0, sizeof(in));
	in[1]=in[7]=0xffff;
	test_fn(in);
	
	printf("--- ALGORITHM STATISTICS ---\n\n");
#define N_TRIES 400000

	clear_stats();
	for (i=0;i<N_TRIES;i++){
		for (ptr=in; ptr<in+12; ptr++)
			*ptr=random();
		ycbcr_square(scratchpad, scratchpad+4, scratchpad+5,
			in, 2);
	}
	printf("ordinary    ");
	print_errors();

	clear_stats();
	for (i=0;i<N_TRIES;i++){
		for (ptr=in; ptr<in+12; ptr++)
			*ptr=random();
		luminaplex_square(scratchpad, scratchpad+4, scratchpad+5,
			in, 2, 1,0);
	}
	printf("Hyperluma 1 ");
	print_errors();

	clear_stats();
	for (i=0;i<N_TRIES;i++){
		for (ptr=in; ptr<in+12; ptr++)
			*ptr=random();
		luminaplex_square(scratchpad, scratchpad+4, scratchpad+5,
			in, 2, 2,0);
	}
	printf("Hyperluma 2 ");
	print_errors();
	clear_stats();
	for (i=0;i<N_TRIES;i++){
		for (ptr=in; ptr<in+12; ptr++)
			*ptr=random();
		luminaplex_square(scratchpad, scratchpad+4, scratchpad+5,
			in, 2, 3,0);
	}
	printf("Luminaplex  ");
	print_errors();

	return 0;
}

void print_usage(void){
	fprintf(stderr,"Usage: %s "
		"_temporal_resample_x_input_width_x_input_height_"
		"x_output_framerate_"
		"x_input_gamma_x_output_width_x_output_height_ "
		"< foo.pix > foo.y4m\n"
		"The input is a concatenation of BRL-CAD .pix files "
		"in the sequence they should appear in the video. "
		"The output is a YUV4MPEG file which is a 4:2:0 "
		"Y'CbCr with simple ASCII headers "
		"that is suitable for transcoding into "
		"a wide variety of video formats, including "
		"Ogg Theora. Input to this program = physical reality"
		"^input_gamma. If you don't know what to write into the"
		" gamma, use 1 if you didn't specify gamma to rt. "
		"If you specified gamma to rt then use 1/that_value. "
		"For example if you specified 2.2 to rt, then use "
		"0.4545 here. It is recommended to specify 2.2 to rt "
		"and then specifying 0.4545 here because that's how you"
		" get the best picture.\n"
		"To be compatible with MPEG, use 1, 5, 10, 12, 15, "
		"24, 25, 30, 50, or 60 framerate. MPEG has other "
		"framerates but these are not integers. \n"
		"If you want PPM output instead of YU4MPEG, prepend "
		"the parameter string with P"
		"(without a space).\n"
		"Temporal resample 1 means no resampling. 2 means"
		" average every 2 consecutive frames into one, "
		"3 average 3 etc."
		" If there are extra frames at the end that wouldn't "
		"make a complete output frame these are discarded.\n"
		"Give a single -table parameter to generate the "
		"hyperluma table to stdout.\n"
		, progname);
}

int main(int argc, char **argv)
{
	float timectr=0;
	unsigned framectr=0;
	unsigned phase;
	float last_hyperluma_switch=0; /* Time of ... */
	double db, rms;
	
	if (argc>=1) progname=(unsigned char *)argv[0];
	if (argc<2){
		print_usage();
		exit(1);
	}
	selftest();

	make_gamma_tables();
	if (!strcmp(argv[1],"-calc_compressed_table")){
		calc_table(1);
		return 0;
	}
	if (!strcmp(argv[1],"-calc_raw_table")){
		calc_table(0);
		return 0;
	}
#ifdef NOTABLE
	fprintf(stderr,"%s: error: hyperluma table is not yet compiled in!\n", progname);
	exit(1);
#else
	uncompress_hyperluma_table();
#endif
	if (!strcmp(argv[1],"-dump_raw_table")){
		fwrite(hyperluma_table, sizeof hyperluma_table, 1, stdout);
		return 0;
	}
	if (!strcmp(argv[1],"-algo-stats")){
		algo_stats();
		return 0;
	}
	if (!strcmp(argv[1],"-testcard")){
		testcard();
		return 0;
	}
	if (!strcmp(argv[1],"-testcard2")){
		desaturate=1;
		testcard();
		return 0;
	}

	if (argv[1][0]=='D'){
		demo=1;
		argv[1]++;
	}
		
	if (argv[1][0]=='P'){
		if (sscanf(argv[1], "P%ux%ux%ux%ux%lfx%ux%u"
				,&temporal_resample
				,&input_w, &input_h, &rate, &input_gamma
				,&output_w
				,&output_h)!=7){
			fprintf(stderr,"%s: Invalid argument format\n", progname);
			print_usage();
			exit(3);
		}
		ppm=1; /* Generate PPM */
	}else{
		if (sscanf(argv[1], "%ux%ux%ux%ux%lfx%ux%u"
				,&temporal_resample
				,&input_w
				,&input_h
				,&rate
				,&input_gamma
				,&output_w
				,&output_h)!=7){
			fprintf(stderr,"%s: Invalid argument format\n", progname);
			print_usage();
			exit(3);
		}
		/* Generate YUV4MPEG2 */
	}
	if (!(input_w&&input_h&&output_w&&output_h&&temporal_resample)){
		fprintf(stderr,"%s: Error: at least one dimension"
				" or temporal resample is zero\n",
				progname);
		exit(7);
	}
	line_buffer=mem_alloc(3*input_w); /* Stays allocated all the time
					     - is small. */

	/* Print the header */
	if (!ppm)
		printf("YUV4MPEG2 W%u H%u F%u:1 Ip A1:1 C420jpeg\n", output_w, output_h, rate);
	fprintf(stderr,"Output format: %s %u fps, temporal resampling: %ux\n",
			ppm?"PPM":"YUV4MPEG2", rate, temporal_resample);
	fprintf(stderr,"Input size: %ux%u, output size: %ux%u\n",
			input_w, input_h, output_w, output_h);

	while(1){
		/* Now d1, d2, d3 are invalid. */

		for (phase=0; phase<temporal_resample; phase++)
			if (load_input(temporal_resample,
						phase)) goto end;
						/* End of input */

		/* Now d1 is valid and d2, d3 invalid. */

#if 0
		/* Logo disabled - annoying */
		if (!ppm&&hyperluma_enable) luminaplex_logo(d1, input_w, input_h);
#endif

		scale_color(d1, input_w, input_h, &d2, output_w, output_h);
		/* Now d1 and d3 are invalid and d2 is allocated */

		clear_stats();
		if (!ppm)dump_luminaplex_frame(d2, output_w, output_h);

		else{
			d2_to_d3();
			/* Now d1 and d2 are invalid and d3 is allocated. */

			if (ppm) dump_ppm_frame(d3);
			else dump_y4m_frame(d3, output_w, output_h);
			mem_free(d3);
			/* Now d1, d2, d3 are invalid. */
		}

		if (ppm){
			fprintf(stderr
				,"Frame %u - %.2f sec. done.\n"
				,framectr, timectr );
		}else{
			rms=rms_err(&db);
			fprintf(stderr
				,"Frame %u"
				" - %.2f sec. done"
				",%.3f usec/pix"
				",%.2f iter/pix"
				",SNR %.2f dB"
				",RMS err %.2f LSB"
				",reuse %.3f%%"
				"\n"
				,framectr, timectr, usec(), iter(), db, rms
				,(double)100*repetitions
					/(repetitions+non_repetitions));
		}

		timectr+=1.0/rate;
		framectr++;
		if (demo
			&& timectr-last_hyperluma_switch>hyperluma_switch_period){
			hyperluma_enable^=1;
			last_hyperluma_switch=timectr;
		}
	}

end:
	/* Now d1, d2, d3 are invalid. */
	free(line_buffer);
	return 0;

}
