Line2D.o: Line2D.cc Line2D.hh Point2D.hh easy_image.hh UsefulFunctions.h
	$(CC) $(CXXFLAGS) -c $< -o $@
