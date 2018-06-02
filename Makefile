default:
	mkdir -p bin
	(cd src; make)
	mv src/maxyee_par .
	mv src/maxyee .
