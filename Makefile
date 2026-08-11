FILE=$(f)
export FILE
default:
	g++ -c -larmadillo -lgsl -lgslcblas -lm -O3 leads_$(FILE).cpp
	mkdir -p ./files
	mkdir -p ./TestCase/Data
	mv *.o ./files
	cp *.h ./files
	cp -ru ./files/ ./TestCase
	cp input.cpp ./TestCase
	cp RunMakefile ./TestCase/Makefile
all:
	g++ -c -larmadillo -lgsl -lgslcblas -lm -O3 leads_*.cpp
	mkdir -p ./files
	mkdir -p ./TestCase/Data
	mv *.o ./files
	cp *.h ./files
	cp -ru ./files/ ./TestCase
	cp input.cpp ./TestCase
	cp RunMakefile ./TestCase/Makefile
clean:
	rm -rf ./files/*
copy:
	cp input.cpp ./TestCase
backup: 
	cp leads*.h leads*.cpp input.cpp /home/Amol/Dropbox/NewProblem/newLEADS/LEADS-2.0
	cp -ru ./Documentation /home/Amol/Dropbox/NewProblem/newLEADS/LEADS-2.0   
	cp Makefile RunMakefile /home/Amol/Dropbox/NewProblem/newLEADS/LEADS-2.0  
	
