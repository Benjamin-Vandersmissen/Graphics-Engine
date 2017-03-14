#include "easy_image.hh"
#include "ini_configuration.hh"
#include "IntroParser.hh"
#include "LSystems.hh"
#include "WireFrameParser.hh"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>



img::EasyImage generate_image(const ini::Configuration &configuration)
{
    std::string type;
    img::EasyImage image;
    type = configuration["General"]["type"].as_string_or_die();
    img::EasyImage image2(10,10);
    image2.draw_line(0,0,9,9,img::Color(100,100,100));

    if (type == "IntroColorRectangle" || type == "IntroBlocks" || type == "IntroLines"){
        IntroParser parser(configuration);
        image = parser.getImage();
    }

    else if (type == "2DLSystem"){
        int size = configuration["General"]["size"].as_int_or_die();
        img::Color BackColor = extractColor(configuration["General"]["backgroundcolor"].as_double_tuple_or_die());
        std::string filename = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        img::Color color = extractColor(configuration["2DLSystem"]["color"].as_double_tuple_or_die());
        bool rainbow = configuration["Extra"]["RainbowLines"].as_bool_or_default(false);
        LParser::LSystem2D system = getLSystem2D(filename);
        image = drawLSystem2D(system, BackColor, color, size, rainbow);
    }

    else if (type == "Wireframe"){
        WireFrameParser parser(configuration);
        image = parser.getImage();
    }

    else if (type == "ZBufferedWireframe"){
        WireFrameParser parser(configuration, true);
        image = parser.getImage();
    }
	return image;
}

int main(int argc, char const* argv[])
{

        int retVal = 0;
        try
        {
                for(int i = 1; i < argc; ++i)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(argv[i]);
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string fileName(argv[i]);
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << argv[i] << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		    //mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}

