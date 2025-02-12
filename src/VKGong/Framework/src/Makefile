CC=gcc
CXX=g++
COMMONDIR=../../common/src/sparse
CFLAGS=-I$(COMMONDIR) -Wno-write-strings -O3 -DUSE_SSE `xml2-config --cflags` -I/home/angels/gsl/include
LDFLAGS=-lpthread `xml2-config --libs` -L/home/angels/gsl/lib -lgsl -lgslcblas -lm
CPPSOURCES=CircularImperfectPlate.cpp Component.cpp GlobalSettings.cpp Input.cpp InputHarmonicCircular.cpp InputHarmonicRectangular.cpp InputNoiseCircular.cpp InputNoiseRectangular.cpp InputStrikeCircular.cpp InputStrikeRectangular.cpp Logger.cpp MathUtil.cpp Output.cpp OutputCircular.cpp OutputRectangular.cpp Parser.cpp PhysicalSystem.cpp PhysicalSystemParser.cpp PhysicalSystemParserCircular.cpp PhysicalSystemParserRectangular.cpp Profiler.cpp RectangularImperfectPlate.cpp ScoreParser.cpp ScoreParserPlate.cpp SettingsManager.cpp WavReader.cpp WavWriter.cpp main.cpp 
CSOURCES=matgen.c
OBJECTS=$(CPPSOURCES:.cpp=.o) $(CSOURCES:.c=.o)
COMMONOBJS=$(COMMONDIR)/csrmatrix.o $(COMMONDIR)/arralloc.o $(COMMONDIR)/pcg.o $(COMMONDIR)/sse.o $(COMMONDIR)/banded.o

all: VKgong

VKgong: $(OBJECTS) $(COMMONOBJS)
	$(CXX) -o VKgong $(OBJECTS) $(COMMONOBJS) $(LDFLAGS)

.cpp.o:
	$(CXX) -c $< -o $@ $(CFLAGS)

.c.o:
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f $(OBJECTS) VKgong

$(COMMONDIR)/csrmatrix.o: $(COMMONDIR)/csrmatrix.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(COMMONDIR)/arralloc.o: $(COMMONDIR)/arralloc.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(COMMONDIR)/pcg.o: $(COMMONDIR)/pcg.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(COMMONDIR)/banded.o: $(COMMONDIR)/banded.c
	$(CC) -c $< -o $@ $(CFLAGS)

$(COMMONDIR)/sse.o: $(COMMONDIR)/sse.c
	$(CC) -c $< -o $@ $(CFLAGS)
