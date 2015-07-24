all: EinsteinExact.out
# Don't create a *.out file until Kranc has finished successfully.
# This way, one can simply re-run "make" when Kranc is interrupted or
# reports an error.
%.out: %.m Metrics/*.m Metrics/*/*.m
	rm -f $*.out
	{ ../../Kranc/Bin/kranc $*.m &&					\
	  (cd .. && git apply m/docfix.patch); } 2>&1 | tee $*.err &&	\
	mv $*.err $*.out
clean:
	rm -f EinsteinExact.out EinsteinExact.err
.PHONY: all clean
