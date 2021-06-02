/* Twibright Pitchotron, a musical transcription analyzer
 * (c) 2007 Karel 'Clock' Kulhavy
 * Twibright Pitchotron is a free software released under the GNU Public License
 * version 3 or later as you wish. */

#define _GNU_SOURCE
/* According to
 * http://www2.informatik.uni-halle.de/lehre/c/libc/libc.texinfo_1.html#SEC12
 * "_GNU_SOURCE If you define this macro, everything is included: ANSI C,
 * POSIX.1, POSIX.2, BSD, SVID, and GNU extensions. In the cases where POSIX.1
 * conflicts with BSD, the POSIX definitions take precedence. "
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define MIN_BRIGHTNESS 0
#define MAX_BRIGHTNESS 255
#define GRATICULE_BRIGHTNESS 64
#define WIDTH 400
#define GRAPH_WIDTH (WIDTH-3*LETTER_WIDTH)
#define LETTER_WIDTH 4 /* Including the space */
#define LETTER_HEIGHT 6 /* Including the space */
#define TOP_SEMITONE 36 /* Above A4 */
#define BOTTOM_SEMITONE -57 /* Above A4 */
#define BANDS (1+TOP_SEMITONE-BOTTOM_SEMITONE)
#define BAND_HEIGHT 6
#if LETTER_HEIGHT != BAND_HEIGHT
#error LETTER_HEIGHT and BAND_HEIGHT must be the same!
#endif
#define HEIGHT (BAND_HEIGHT*BANDS)
#define FRAME_RATE 12 /* fps */
#define ANALYSIS_DYNAMIC_RANGE 22 /* dB */
#define NOVERT 7 /* How many overtones to take into account. */

#define VISIBLE_TIME 4
#define ALC_DECAY 5 /* decibels per second. Must be nonnegative. */
/* How many seconds of audio are on the screen at a time */

float middle_a=440;

struct osc{
	double x,y;
	double rot_x, rot_y;
	double decay;
	double energy_scale;
};

int weighted_overtones=1;
unsigned char array[GRAPH_WIDTH*(BANDS+1)]; /* The lowest is for the graticule */
int xpos; /* First unwritten. */
float results[BANDS]; /* Cumulative energy */
float summed_results[BANDS];
float export_results[BANDS];
struct osc oscs[BANDS];
double pixel_period, sample_period; /* Periods in seconds */
unsigned pixel_divisor, idle_frames_countdown;
unsigned long idle_samples; /* idle_frames_countdown expressed in samples 
			       instead of frames */
unsigned char scroll; /* 1=scroll, 0=leave picture in place */
unsigned char removal=1; /* Overtone removal */
unsigned char output_summed=1; /* 1=outputs summed, 0=outputs the original one */
float removal_dummy=0; /* 0 means remove before logarithmization, -HUGE_VAL
			  after */
float relative_bandwidth=1;
int n_channels;
float sample_rate;
int verbose;
double bpm; /* How fast to print the graticule. If zero, no graticule is
	       printed. */
double pixel_graticules; /* How many graticules one pixel represents - this
			    number is typically much less than one */

unsigned char label[3*3*LETTER_WIDTH*BAND_HEIGHT*BANDS]; /* Only for one
							 octave */
float display_dynamic_range=20; /* dB */
float zero_guard=90; /* If automatic ranging and overtone removal are both on,
			then the overtone removal produces zeroes which map
			to -Inf and that kills the dynamic range - then
			everything is black. With zero_guard, decibel
			values less than max=zero_guard are ignored and
			not taken into account for calculation of the
			dynamic range at all. */

/* Begins with C - red */
static unsigned short colours[13*3]={
	257,  0,  0,
	  0,160,257,
	257,257,  0,
	180,  0,257,
	  0,257,  0,
	257,  0,115, 
	  0,257,257,
	257,140,  0,
	  0,  0,257,
	190,257,  0,
	257,  0,257, 
	  0,257,160,
	257,257,257 /* This is for the graticule */
	/* 257 is here instead of 256 so that we don't have to add
	 * 0x80 before the shift, hence speeding it up a bit */
};

const unsigned char letters[]={
	/* Space */
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,

	/* C */
	0,1,1,0,
	1,0,0,0,
	1,0,0,0,
	1,0,0,0,
	0,1,1,0,
	0,0,0,0,

	/* D */
	1,1,0,0,
	1,0,1,0,
	1,0,1,0,
	1,0,1,0,
	1,1,0,0,
	0,0,0,0,

	/* E */
	1,1,1,0,
	1,0,0,0,
	1,1,1,0,
	1,0,0,0,
	1,1,1,0,
	0,0,0,0,

	/* F */
	1,1,1,0,
	1,0,0,0,
	1,1,0,0,
	1,0,0,0,
	1,0,0,0,
	0,0,0,0,

	/* G */
	0,1,1,0,
	1,0,0,0,
	1,0,0,0,
	1,0,1,0,
	0,1,1,0,
	0,0,0,0,

	/* A */
	0,1,0,0,
	1,0,1,0,
	1,1,1,0,
	1,0,1,0,
	1,0,1,0,
	0,0,0,0,

	/* B */
	1,1,0,0,
	1,0,1,0,
	1,1,0,0,
	1,0,1,0,
	1,1,0,0,
	0,0,0,0,

	/* # */
	1,0,1,0,
	1,1,1,0,
	1,0,1,0,
	1,1,1,0,
	1,0,1,0,
	0,0,0,0,

	/* 0 */
	0,1,0,0,
	1,0,1,0,
	1,0,1,0,
	1,0,1,0,
	0,1,0,0,
	0,0,0,0,

	/* 1 */
	0,1,0,0,
	1,1,0,0,
	0,1,0,0,
	0,1,0,0,
	0,1,0,0,
	0,0,0,0,

	/* 2 */
	0,1,0,0,
	1,0,1,0,
	0,0,1,0,
	0,1,0,0,
	1,1,1,0,
	0,0,0,0,

	/* 3 */
	1,1,0,0,
	0,0,1,0,
	0,1,0,0,
	0,0,1,0,
	1,1,0,0,
	0,0,0,0,

	/* 4 */
	0,0,1,0,
	0,1,1,0,
	1,0,1,0,
	1,1,1,0,
	0,0,1,0,
	0,0,0,0,

	/* 5 */
	1,1,1,0,
	1,0,0,0,
	1,1,0,0,
	0,0,1,0,
	1,1,0,0,
	0,0,0,0,

	/* 6 */
	0,1,1,0,
	1,0,0,0,
	1,1,0,0,
	1,0,1,0,
	0,1,0,0,
	0,0,0,0,

	/* 7 */
	1,1,1,0,
	0,0,1,0,
	0,1,0,0,
	0,1,0,0,
	0,1,0,0,
	0,0,0,0,

	/* 8 */
	0,1,0,0,
	1,0,1,0,
	0,1,0,0,
	1,0,1,0,
	0,1,0,0,
	0,0,0,0,

	/* 9 */
	0,1,0,0,
	1,0,1,0,
	0,1,1,0,
	0,0,1,0,
	1,1,0,0,
	0,0,0,0,

};

int overtones[NOVERT]; /* Semitone offsets for 1th, 2th, 3th etc. harmonic */

static void array_write(const unsigned char *ary, unsigned start, unsigned len, int col){
	unsigned short *colptr=colours+col*3;

	ary+=start;
	for (;len;len--){
		putchar((*ary*colptr[0])>>8);
		putchar((*ary*colptr[1])>>8);
		putchar((*ary*colptr[2])>>8);
		ary++;
	}
}

static void array_write_wrap(const unsigned char *ary, unsigned start, unsigned len, int col){
	while(len){
		unsigned consume=len;
		if (start+consume>GRAPH_WIDTH) consume=GRAPH_WIDTH-start;
		array_write(ary, start, consume, col);
		start+=consume;
		if (start>=GRAPH_WIDTH) start=0;
		len-=consume;
	}
}

/* semitone can be 0...12 - the semitone to display
 * curpos is the cursor (vertical bar) position. */
/* curpos 0...GRAPH_WIDTH-1 */
static void dump_microline(unsigned char *src, int semitone, int curpos)
{
	static const unsigned char full=MAX_BRIGHTNESS;

	if (scroll){
		int w=GRAPH_WIDTH>>1;
		array_write_wrap(src, xpos, w-1, semitone);
		array_write(&full, 0, 1, semitone);
		curpos=xpos+w;
		if (curpos>=GRAPH_WIDTH) curpos-=GRAPH_WIDTH;
		array_write_wrap(src, curpos, GRAPH_WIDTH-w, semitone);
	}else{
		array_write(src, 0, curpos, semitone);
		array_write(&full, 0, 1, semitone);
		if (curpos<GRAPH_WIDTH-1) array_write(src,curpos+1, GRAPH_WIDTH-curpos-1,semitone);
	}
}

static void dump_frame(void)
{
	int i, band;
	int curpos;
	/* cursor position, half cycle apart from xpos */
	unsigned char *labelptr;

	curpos=xpos+(GRAPH_WIDTH>>1); 
	if (curpos>=GRAPH_WIDTH) curpos-=GRAPH_WIDTH;
	labelptr=label;
	for (band=0;band<BANDS;band++){
		int semitone=(BOTTOM_SEMITONE+band-3)%12;
		if (semitone<0) semitone+=12;
		for (i=0;i<BAND_HEIGHT;i++){
			int display=(!i||i==BAND_HEIGHT-1)?12:semitone;
			fwrite(labelptr, 3*LETTER_WIDTH, 3, stdout);
			labelptr+=3*LETTER_WIDTH*3;
			dump_microline(array+(display==12?BANDS:band)*GRAPH_WIDTH
				, display, curpos);
		}
	}
}

int find_strongest_summed_results(void)
{
	int i, maxpos;
	float max;

	maxpos=0;
	max=*summed_results;

	for (i=1;i<sizeof summed_results/sizeof *summed_results; i++)
		if (summed_results[i]>max){
			max=summed_results[i];
			maxpos=i;
		}
	return maxpos;
}

void remove_halftone(int halftone){
	int overtone;

	summed_results[halftone]=removal_dummy;
	for (overtone=0;overtone<NOVERT;overtone++){
		int target=halftone-overtones[overtone];
		if (target<0) break;
		if (summed_results[target])summed_results[target]
			-=weighted_overtones?
				results[halftone]*(float)(NOVERT-overtone)/NOVERT
				:results[halftone];
		/* Can happen that it's already zeroed out. */
	}

}

void remove_overtones(void){
	int halftone, overtone;

	for (halftone=0;halftone<BANDS;halftone++){
		float sum;
		sum=results[halftone];

		for (overtone=0;overtone<NOVERT;overtone++){
			int target=halftone+overtones[overtone];
			if (target>=BANDS) break;
			sum+=weighted_overtones?
				results[target]*(float)(NOVERT-overtone)/NOVERT
				:results[target];
			/* If NOVERT is 4, the weights will be 1, 0.75, 0.5,
			 * 0.25 */
		}
		/* Sum the energy together */
		summed_results[halftone]=sum;

	}
	/* Now results2 have their fundamental and overtone energies summed. */
	memset(export_results, 0, sizeof results);

	while(summed_results[halftone=find_strongest_summed_results()]!=removal_dummy){
		export_results[halftone]=output_summed?
			summed_results[halftone]:results[halftone];
		remove_halftone(halftone);
		for (overtone=0;overtone<NOVERT;overtone++){
			int target=halftone+overtones[overtone];
			if (target>=BANDS) break;
			remove_halftone(target);
		}
	}
	memcpy(results, export_results, sizeof export_results);
}

static void get_results(void)
{
	float *r;
	static float max=-HUGE_VAL; /* Maximum decibels */
	float min;

	/* Now the results mean energy */
	if (removal_dummy==0&&removal){
		/* Perform removal on energy */
		remove_overtones();
	}

	for (r=results;r<results+sizeof results/sizeof *results; r++)
		*r=10*log10(*r); /* Logarithmize */

	/* Now the results mean decibels, the higher the stronger */
	if (removal_dummy==-HUGE_VAL&&removal){
		/* Perform removal on decibel (cepstrum) */
		remove_overtones();
	}
		
	/* Update the maximum expressed always in decibels */
	max-=(float)ALC_DECAY/pixel_period;
	for (r=results;r<results+sizeof results/sizeof *results; r++)
		if (*r>max) max=*r;

	if (display_dynamic_range){
		/* Pre-selected dynamic range */
		min=max-display_dynamic_range;
	}else{
		/* Dynamic range automatically calculated to fit minimum and
		 * maximum current value */

		/* Find the minimum */
		min=max;
		for (r=results;r<results+sizeof(results)/sizeof(*results);
				r++)
			if (*r>=max-zero_guard&&*r<min) min=*r;
			/* Ignore values that are zero amplitude (-Inf dB)
			 * or extremely small values. */
	}
	for (r=results;r<results+sizeof results/sizeof *results; r++){
		*r=floor((*r-min)*(MAX_BRIGHTNESS-MIN_BRIGHTNESS)/
				(display_dynamic_range?display_dynamic_range:
				 (max-min))
				+0.5);
		if (*r<0) *r=0;
		*r+=MIN_BRIGHTNESS;
	}
}

static void init_oscs(void)
{
	int osc;
	float neighbour_response; /* POWER response of the lower
				     neighbouting semitone */
	float rc_relative_freq;
	neighbour_response=pow(10,-ANALYSIS_DYNAMIC_RANGE/10);
	if (verbose){
		fprintf(stderr,"Neighbour response %G dB means amplitude "
			"response %G., ", (double)-ANALYSIS_DYNAMIC_RANGE
			, (double)neighbour_response);
	}
	rc_relative_freq
		=sqrt(1/neighbour_response-1);
	if (verbose){
		fprintf(stderr,"that's relative frequency %G of the "
			"lowpass filter cutoff frequency.\n", rc_relative_freq);
	}

	/* osc=0 is BOTTOM_SEMITONE etc. */
	for (osc=0;osc<sizeof oscs/sizeof *oscs;osc++){
		float freq;
		float angle;
		float half_bandwidth;
		float sine_response; /* What approx. amplitude (not energy) the
					oscillator will have when a sine wave of
					the center frequency is present */
		float decay_time; /* After this time the amplitude is e times
				     smaller and energy e^2 times smaller */

		freq=middle_a*pow(2,(osc+BOTTOM_SEMITONE)/12.0);
		if (verbose)
			fprintf(stderr," Semitone %d has frequency %G Hz, ",
				osc, freq);

		/* Lowpass filter bandwidth */
		half_bandwidth=freq*(1-pow(2,-1.0/12)); 

		half_bandwidth*=relative_bandwidth;

		/* Now if relative_bandwidth is 1 then half_bandwidth is that 
		 * the -3dB (half energy) point
		 * lies on the lower neighbouring semitone. */
		if (verbose)
			fprintf(stderr," distance to neighbouring lower "
				"semitone %G Hz, ", half_bandwidth);

		half_bandwidth/=rc_relative_freq;
		decay_time=1/(half_bandwidth*2*M_PI);
		if (verbose)
			fprintf(stderr," half bandwidth %G Hz, amplitude "
				"exponential decay time %G ms and "
				"decays through the requested dynamic "
				"range %G dB in %G ms\n"
				,half_bandwidth
				,1000*decay_time
				,(double)ANALYSIS_DYNAMIC_RANGE
				,1000*(decay_time*ANALYSIS_DYNAMIC_RANGE
					/(10*log10(exp(2)))));
		angle=2*M_PI*freq/sample_rate;
		oscs[osc].x=0;
		oscs[osc].y=0;
		oscs[osc].rot_x=cos(angle);
		oscs[osc].rot_y=sin(angle);
		oscs[osc].decay=exp(-1/(sample_rate*decay_time));

		/* If the frequency were zero, the oscillator would behave
		 * as a lowpass RC filter and the sum of convolving a constant
		 * function y=1 and the RC filter would be given by a
		 * harmonic series: */
		sine_response=1/(1-oscs[osc].decay);

		/* But the damped oscillator behaves like a RC lowpass
		 * multiplied by a sine wave of the filter zero crossing
		 * frequency.
		 *
		 * To find out the response to the zero crossing frequency,
		 * we use the theorem about frequency of product: it's
		 * a convolution of the spectra. And the DC coefficient of the
		 * RC lowpass is the sum of harmonic series we have in
		 * sine_response. */

		/* If two sine waves multiply, they result in ripples that
		 * have only 0.5 average height. */
		sine_response/=2;
		oscs[osc].energy_scale=1/(sine_response*sine_response);
	}
}

/* Returns oscillator energy after the step. */
static float step_osc(struct osc *o, float input)
{
	double x,y;

	x=o->x*o->rot_x-o->y*o->rot_y;
	y=o->y*o->rot_x+o->x*o->rot_y;
	x*=o->decay;
	y*=o->decay;
	x+=input;

	o->x=x;
	o->y=y;
	return (x*x+y*y)*o->energy_scale;
}

static void step_oscs(float sample)
{
	struct osc *o;
	for (o=oscs; o<oscs+sizeof oscs/sizeof *oscs; o++)
		results[o-oscs]+=step_osc(o, sample);
}

static void next_column(void)
{
	int band;
	static double graticule_remainder; /* 0<=graticule_remainder<1 */

	get_results();
	for (band=0;band<BANDS;band++){
		array[band*GRAPH_WIDTH+xpos]=results[band];
	}
	graticule_remainder+=pixel_graticules;
	if (graticule_remainder>=1){
		graticule_remainder-=1;
		array[BANDS*GRAPH_WIDTH+xpos]=GRATICULE_BRIGHTNESS;
	}else array[BANDS*GRAPH_WIDTH+xpos]=0;
	xpos++;
	if (xpos>=GRAPH_WIDTH) xpos=0;
	memset(results, 0, sizeof(results)); 
	/* Prepare for the next column summation */
}

/* Signed is important here - the WAV data are signed! */
void process_input_sample(signed char *samples)
{
	float sum;

	sum=(unsigned char)samples[0]+((signed char)samples[1]<<8);
	if (n_channels==2)
		sum+=(unsigned char)samples[2]+((signed char)samples[3]<<8);

	step_oscs(sum/n_channels/32767.5);
}

void process_pixels(void)
{
	static double time_remainder;
	static unsigned pixel_remainder;
	static double time;

	time_remainder+=sample_period;
	if (time_remainder>=pixel_period){
		/* Time for writing a pixel */
		time_remainder-=pixel_period;
		next_column();
		pixel_remainder++;
		if (pixel_remainder>=pixel_divisor){
			/* Time for dumping a frame */
			pixel_remainder=0;
			if (idle_frames_countdown) idle_frames_countdown--;
			else{
				dump_frame();
				time+=pixel_period*pixel_divisor;
				fprintf(stderr,"Frame %lG sec.\n", time);
			}
		}
	}
}

void init_harmonics(void)
{
	int h;

	for (h=0;h<NOVERT;h++){
		overtones[h]=floor(12*log(h+2)/log(2)+0.5);
		fprintf(stderr,"%dth harmonic %d semitones up\n", h+2
				,overtones[h]);
	}
}

/* 0-7 are letters C D E F G A B, 8 is #, 9-18 are digits 0-9 */
static void print_char(unsigned char *dest, char c, unsigned short *colour)
{
	const unsigned char*src
		=letters+LETTER_WIDTH*((c+1)*LETTER_HEIGHT-1);
	int x,y;

	for (y=0;y<LETTER_HEIGHT;y++){
		for (x=0;x<LETTER_WIDTH;x++){
			if (*src++){
				dest[0]=(colour[0]*255)>>8;
				dest[1]=(colour[1]*255)>>8;
				dest[2]=(colour[2]*255)>>8;
			};
			dest+=3;
		}
		src-=LETTER_WIDTH<<1;
		dest+=3*LETTER_WIDTH*2;
	}
}

static void print(unsigned char *dest, unsigned char *text, unsigned short
	*colour)
{
	int c;

	for (c=3;c;c--){
		print_char(dest, *text, colour);
		dest+=3*LETTER_WIDTH;
		text++;
	}
}

void make_label(void)
{
	unsigned char *labelptr;
 	int band;

	labelptr=label;
	for (band=0;band<BANDS;band++){
		unsigned char txt[3];
		int semitone=(BOTTOM_SEMITONE+band-3)%12;
		int octave=(72+BOTTOM_SEMITONE+band-3)/12-1;
		if (semitone<0) semitone+=12;
		txt[0]=semitone;
		if (semitone>=5) txt[0]++;
		txt[1]=txt[0]&1?8:0;
		txt[0]>>=1;
		txt[0]++;
		if (octave>=9||octave<0) octave=-1;
		txt[2]=octave+9;
		print(labelptr, txt, colours+3*semitone);
		labelptr+=BAND_HEIGHT*3*3*LETTER_WIDTH;
	}
}

void init(void)
{
	pixel_period=(double)VISIBLE_TIME/GRAPH_WIDTH;
	sample_period=(double)1/sample_rate;
	memset(array, MIN_BRIGHTNESS, sizeof array);
	{
		double frame_period;

		frame_period=(double)1/FRAME_RATE;
		pixel_divisor=floor(frame_period/pixel_period+0.5);
		pixel_period=frame_period/pixel_divisor;
		idle_frames_countdown=floor(GRAPH_WIDTH*pixel_period/2
			/frame_period+0.5);
		idle_samples=floor(GRAPH_WIDTH*pixel_period/2*sample_rate+0.5);

	}
	init_oscs();
	init_harmonics();
	pixel_graticules=pixel_period*bpm/60;
	fprintf(stderr,"pixel period=%G, pixel graticules=%G\n"
			, pixel_period, pixel_graticules);
	make_label();
}

void parse_args(int argc, char **argv)
{
next_opt:
	switch(getopt(argc, argv, "hsolb:r:a:d:")){
		case 'h':
			fprintf(stderr,"Twibright Pitchotron, a musical "
				"transcription pitch analyzer.\n"
				"Usage: pitchotron-bin <arguments> < music.wav\n"
				"\n"
				"-h print help\n"
				"-s scroll display, in place otherwise.\n"
				"-o do not remove overtones\n"
				"-l remove after logarithmization, "
					"otherwise remove before\n"
				"-b relative_bandwidth Set relative "
					"bandwidth. 1 means the lower "
					"neighbouring halftone excites the "
					"halftone oscillator at -3dB. Defaults"
					" to 1.\n"
				"-r rhythm_in_bpm sets the vertical graticule"
					"\n" 
				"-d set the display dynamic range in dB. "
					" Default 20dB. Must not be negative. "
					"If zero, the range is automatically "
					"fitted at every moment.\n"
				"-a freq_in_Hz Set middle A in Hz.\n"
				);
					exit(1);
			break;
		case 's': scroll=1; break;
		case 'o': removal=0; break;
		case 'l': removal_dummy=-HUGE_VAL; break;
		case 'b': relative_bandwidth=atof(optarg); break;
		case 'r': bpm=atof(optarg); break;
		case 'a': middle_a=atof(optarg); break;
		case 'd': display_dynamic_range=atof(optarg);
			  if (display_dynamic_range<0){
				  fprintf(stderr,"Display dynamic range "
					"cannot be negative!\n");
				  exit(1);
			  }
			  break;
		case -1:
			return;
	}
	goto next_opt;
}

/* Wastes "len" chars from the stdin */
void waste(unsigned len){
	for (;len;len--) getchar();
}

void read_wav_header(void){
	static unsigned char buf[44];

	if (fread(buf, 44, 1, stdin)!=1){
		fprintf(stderr,"pitchotron: error: cannot load WAV header.\n");
		exit(1);
	}
	if (memcmp(buf,"RIFF",4)){
not_a_wav:
		fprintf(stderr,"pitchotron error: input not a PCM wav file.\n");
		exit(1);
	}
	if (memcmp(buf+8,"WAVEfmt \x10\0\0\0",12)) goto not_a_wav;
	n_channels=buf[22];
	if (buf[23]||n_channels>2||n_channels<1){
		fprintf(stderr,"pitchotron: error: only 1-2 channels in wav"
				" file supported\n");
		exit(1);
	}
	sample_rate=buf[24]+((unsigned)buf[25]<<8)+((unsigned long)buf[26]<<16)
		+((unsigned long)buf[27]<<24);
	if (buf[34]!=16||buf[35]!=0){
		fprintf(stderr,"pitchotron: error: only 16 bits per sample WAVs"
			" supported.\n");
		exit(1);
	}

	fprintf(stderr,"pitchotron input: sample rate "
		"%G kHz, %u channels\n", sample_rate/1000, n_channels);
}

int main(int argc, char**argv)
{
	static signed char samples[4];

	parse_args(argc, argv);
	read_wav_header();
	init();

	while (fread(samples,2,n_channels,stdin)==n_channels){
		process_input_sample(samples);
		process_pixels();
	}
	for (;idle_samples;idle_samples--){
		step_oscs(0);
		process_pixels();

	}
	return 0;
}
