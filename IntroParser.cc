//
// Created by Benjamin on 14/02/2017.
//
#include "IntroParser.hh"

void IntroParser::parseColorRectangle(const ini::Configuration &configuration) {
    int w = configuration["ImageProperties"]["width"].as_int_or_die();
    int h = configuration["ImageProperties"]["height"].as_int_or_die();
    bool scale = configuration["ColorProperties"]["scaled"].as_bool_or_default(false);
    img::EasyImage image(w,h);

    if (!scale) {
        for (unsigned int i = 0; i < w; i++) {
            for (unsigned int j = 0; j < h; j++) {
                image(i, j).red = i % 256;
                image(i, j).green = j % 256;
                image(i, j).blue = (i + j) % 256;
            }
        }
    }else{
        for(float i = 0; i < w; i++){
            for (float j = 0; j < h; j++){
                image(i,j).red = (int)255*(i/w);
                image(i,j).green = (int)255*(j/h);
                image(i,j).blue = (int)(255*(i/w + j/h)) % 256;
            }
        }
    }
    this->image = image;
}

void IntroParser::parseBlocks(const ini::Configuration &configuration) {
    int Wi = configuration["ImageProperties"]["width"].as_int_or_die();
    int Hi = configuration["ImageProperties"]["height"].as_int_or_die();
    int nrXBlocks = configuration["BlockProperties"]["nrXBlocks"].as_int_or_die();
    int nrYBlocks = configuration["BlockProperties"]["nrYBlocks"].as_int_or_die();
    std::vector<int> ColorW = extractColor(configuration["BlockProperties"]["colorWhite"].as_double_tuple_or_die());
    std::vector<int> ColorB = extractColor(configuration["BlockProperties"]["colorBlack"].as_double_tuple_or_die());
    bool invert = configuration["BlockProperties"]["invertColors"].as_bool_or_default(false);
    img::EasyImage image(Wi,Hi);

    float Wb = Wi / nrXBlocks;
    float Hb = Hi / nrYBlocks;
    for(int i = 0; i < Wi; i++){
        for(int j = 0; j < Hi; j++){
            int Bx = i / Wb;
            int By = j / Hb;
            if (((Bx+By) % 2 == 0 && !invert) ||((Bx+By) % 2 == 1 && invert)){
                image(i,j) = img::Color(ColorB[0], ColorB[1], ColorB[2]);
            }
            else{
                image(i,j) = img::Color(ColorW[0], ColorW[1], ColorW[2]);
            }
        }
    }
    this->image = image;
}

void IntroParser::parseLines(const ini::Configuration &configuration) {
    int w = configuration["ImageProperties"]["width"].as_int_or_die();
    int h = configuration["ImageProperties"]["height"].as_int_or_die();
    std::string figure = configuration["LineProperties"]["figure"].as_string_or_die();
    std::vector<int> ColorBG = extractColor(configuration["LineProperties"]["bgColor"].as_double_tuple_or_die());
    std::vector<int> ColorLine = extractColor(configuration["LineProperties"]["lineColor"].as_double_tuple_or_die());
    int nrLines = configuration["LineProperties"]["nrLines"].as_int_or_die();

    img::EasyImage image(w,h);
    this->image = image;
    image.clear(img::Color(ColorBG[0], ColorBG[1], ColorBG[2]));

    if (figure == "Eye") {
        this->Eye(w, h, ColorLine, nrLines);
    }
    else if (figure == "QuarterCircle") {
        this->QuarterCircle(w, h, ColorLine, nrLines);
    }
    else if (figure == "Diamond"){
        this->Diamond(w, h, ColorLine, nrLines);
    }
}

IntroParser::IntroParser(const ini::Configuration &configuration) {
    std::string type = configuration["General"]["type"].as_string_or_die();
    if (type == "IntroColorRectangle"){
        this->parseColorRectangle(configuration);
    }
    else if (type == "IntroBlocks"){
        this->parseBlocks(configuration);
    }
    else if (type == "IntroLines"){
        this->parseLines(configuration);
    }
}

void IntroParser::QuarterCircle(unsigned int w, unsigned int h, std::vector<int> ColorLine, unsigned int nrLines,
                                unsigned int quadrant, unsigned int x, unsigned int y) {
    int x1 = x; //positie 1 is op CONSTANTE x=0 of x = IMG_SIZEX
    int y1 = y;
    int x2 = x; //positie 2 is op CONSTANTE y=0 of y = IMG_SIZEY
    int y2 = y;
    float d1 = 1; //delta(y) waarmee positie 1 toeneemt per loop
    float d2 = 1; //delta(x) waarmee positie 2 toeneemt per loop

    img::Color colLine(ColorLine[0],ColorLine[1],ColorLine[2]);
    switch(quadrant){
        case 1:
            // y1 = y ~ d1 = 1
            x1 += w-1;
            x2 += w-1;
            y2 += h-1;
            d2 = -1;
            break;

        case 2:
            //x1 = x ~ x2 = x ~ y1 = y ~ d1 = 1 ~ d2 = 1
            y2 += h - 1;
            break;
        case 3:
            //x1 = x ~ x2 = x ~ y2 = y ~ d2 = 1
            y1 += h - 1;
            d1 = -1;
            break;

        case 4:
            //y2 = y
            x1 += w-1;
            y1 += h-1;
            x2 += w-1;
            d1 = -1;
            d2 = -1;
            break;
    }

    if (nrLines != 1) {
        d1 *= (h / (nrLines - 1));
        d2 *= (w / (nrLines - 1));
    }

    for(int i = 0; i < nrLines; i++){
        this->image.draw_line(x1, y1+i*d1 , x2+i*d2, y2, colLine);
    }
    this->image.draw_line(x1, d1 > 0 ? y1+h-1 : y, d2 > 0 ? x2+w-1 : x, y2, colLine);
}

void
IntroParser::Eye(unsigned int w, unsigned int h, std::vector<int> ColorLine, unsigned int nrLines) {
    QuarterCircle(w, h, ColorLine, nrLines, 2);
    QuarterCircle(w, h, ColorLine, nrLines, 4);
}

void IntroParser::Diamond(unsigned int w, unsigned int h, std::vector<int>ColorLine, unsigned int nrLines){
    QuarterCircle(w/2,h/2, ColorLine,nrLines, 1);
    QuarterCircle(w/2,h/2, ColorLine,nrLines, 2, w/2);
    QuarterCircle(w/2,h/2, ColorLine,nrLines, 3, w/2, h/2);
    QuarterCircle(w/2,h/2, ColorLine,nrLines, 4, 0, h/2);

}

const img::EasyImage &IntroParser::getImage() const {
    return image;
}
