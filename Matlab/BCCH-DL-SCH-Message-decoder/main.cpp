#include <stdio.h>
#include <sys/types.h>
#include <BCCH-DL-SCH-Message.h>
/* Rectangle ASN.1 type */
int main(int ac, char **av) {
char buf[1024];
/* Temporary buffer
*/
asn_dec_rval_t rval; /* Decoder return value */
Rectangle_t *rectangle = 0; /* Type to decode. Note this 01 ! */
FILE *fp;
size_t size;
char *filename;
/* Input file handler
/* Number of bytes read
/* Input file name */
*/
*/
/* Require a single filename argument */
if(ac != 2) {
fprintf(stderr, ”Usage: %s <file.ber>\n”, av[0]);
exit(1);
} else {
filename = av[1];
}
/* Open input file as read-only binary */
fp = fopen(filename, ”rb”);
if(!fp) {
perror(filename);
exit(1);
}
/* Read up to the buffer size */
size = fread(buf, 1, sizeof(buf), fp);
fclose(fp);
if(!size) {
fprintf(stderr, ”%s: Empty or broken\n”, filename);
exit(1);
}
/* Decode the input buffer as Rectangle type */
rval = ber_decode(0, &asn_DEF_Rectangle, (void **)&rectangle, buf, size);
if(rval.code != RC_OK) {
fprintf(stderr, ”%s: Broken Rectangle encoding at byte %ld\n”, filename,
(long)rval.consumed);
exit(1);
}
/* Print the decoded Rectangle type as XML */
xer_fprint(stdout, &asn_DEF_Rectangle, rectangle);
return 0; /* Decoding finished successfully */
}
