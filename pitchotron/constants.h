#define TV_GAMMA 2.2


/* Y709 = 0.2125 R + 0.7154 G + 0.0721 B
 * "Contemporary CRT phosphors are standardized in Rec. 709 "
 * (http://www.poynton.com/PDFs/ColorFAQ.pdf) 
 * Conversion from RGB to luminance (Y) 
 * They must make exactly 65536 together. Must be UL because sometimes they
 * sum to almost the whole 32-bit range. */
#define LUMI_R_16 13927UL /*  0.2125, inced to get 65536 total */
#define LUMI_G_16 46884UL /*  0.7154 */
#define LUMI_B_16  4725UL /*  0.0721 */


/* Converstion from 0-255 R'G'B' to 16-235,16-240 Y'CbCr According to BT.601
 * (non-HDTV) */
/* Must sum to 65536 */
#define LUMA_R_16 16829L /* 0.2567 */
#define LUMA_G_16 33039L /* 0.5041 */
#define LUMA_B_16  6416L /* 0.0979 */

/* Must sum to 0 */
#define CB_R_16  -9714L /* -0.148 */
#define CB_G_16 -19070L /* -0.291 */
#define CB_B_16  28784L /*  0.4392*/

/* Must sum to 0 */
#define CR_R_16  28784L /*  0.4392 */
#define CR_G_16 -24103L /* -0.367 */
#define CR_B_16  -4681L /* -0.071 */


#define LUMA_16   76284L /*  1.164 */
#define R_CR_16  104595L /*  1.596 */
#define G_CR_16  -53281L /* -0.813 */
#define G_CB_16  -25690L /*  0.392 */
#define B_CB_16  132186L /*  2.017 */
