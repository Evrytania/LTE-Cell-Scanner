36331-850.asn is genrated as doc/How_to_use_asn1c_lameditor.txt

follow asn1c-usage.pdf, you need test.asn1 and main.c
chapter 3.1 and 3.2


--------------for LTE asn1c decode------------------------------

a ref:
http://blog.csdn.net/peng_yw/article/details/22437251

cd LTE-BCCH-DL-SCH-decode
asn1c  -S /usr/local/share/asn1c -fcompound-names -fskeletons-copy -gen-PER -pdu=auto 36331-850.asn

converter-sample.c
add 
#define PDU BCCH_BCH_Message
#define ASN_PDU_COLLECTION
after 
#include <asn_internal.h>

per_opentype.c

add 
padding = padding % 8;
after ASN_DEBUG("Too large padding %d in open type", (int)padding);
and comment out following:
_ASN_DECODE_FAILED;

make -f Makefile.am.sample

exe:
progname
