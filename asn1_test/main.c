#include <stdio.h>
#include <sys/types.h>
#include <Rectangle.h>
int main(int ac, char **av) {
  char buf[1024];
  asn_dec_rval_t rval;
  Rectangle_t *rectangle = 0;
  FILE *fp;
  size_t size;
  char *filename;

  if(ac != 2){
    fprintf(stderr, "Usage: %s <file.ber>\n", av[0]);
    exit(1);
  } else {
    filename = av[1];
  }

  fp = fopen(filename, "rb");
  if(!fp) {
    perror(filename);
    exit(1);
  }

  size = fread(buf, 1, sizeof(buf), fp);
  fclose(fp);
  if(!size) {
    fprintf(stderr, "%s: Empty or broken\n", filename);
    exit(1);
  }

  rval = ber_decode(0, &asn_DEF_Rectangle, (void **)&rectangle, buf, size);
  if(rval.code != RC_OK) {
    fprintf(stderr, "%s: Broken Rectangle encoding at byte %ld\n", filename,
    (long)rval.consumed);
    exit(1);
  }

  xer_fprint(stdout, &asn_DEF_Rectangle, rectangle);

  return 0;
}

