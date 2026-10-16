CC=gcc
CFLAGS = -W -Wall -I$(ASCDS_INSTALL)/include
LDFLAGS = -L$(ASCDS_INSTALL)/ots/lib -L$(ASCDS_LIB) -lardlib -lascdm -lcaldb4 -lNewHdr2 -lcxcparam -lregion -lreadline -lhistory -lstk -lcfitsio -lm

ardlib_qe.linux : ardlib_qe.c
	$(CC) $(CFLAGS) -v -o $(@) $< $(LDFLAGS)

ardlib_qe.darwin : %.c
	$(CC) $(CFLAGS) -o $(@).darwin $< $(LDFLAGS)

% : %.c
	$(CC) $(CFLAGS) -o $(@) $< $(LDFLAGS)

#%.o : %.c
#	$(CC) -c $(CFLAGS) $< -o $@
