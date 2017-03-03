//
// Created by Benjamin on 14/02/2017.
//
#include "Intro.hh"


img::EasyImage ColorRectangle(unsigned int w, unsigned int h, bool scale){
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
    return image;
}


img::EasyImage Blocks(unsigned int Wi, unsigned int Hi, unsigned int nrXBlocks, unsigned int nrYBlocks, std::vector<int> ColorW, std::vector<int> ColorB, bool invert){
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
    return image;
}

img::EasyImage QuarterCircle(unsigned int w, unsigned int h, std::string figure, std::vector<int> ColorLine, unsigned int nrLines, img::EasyImage& image, unsigned int quadrant = 2, unsigned int x = 0, unsigned int y = 0){
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
        image.draw_line(x1, y1+i*d1 , x2+i*d2, y2, colLine);
    }
    image.draw_line(x1, d1 > 0 ? y1+h-1 : y, d2 > 0 ? x2+w-1 : x, y2, colLine);
    return image;
}

img::EasyImage Eye(unsigned int w, unsigned int h, std::string figure, std::vector<int> ColorBG, std::vector<int>ColorLine, unsigned int nrLines){
    img::EasyImage image(w,h);
    image.clear(img::Color(ColorBG[0], ColorBG[1], ColorBG[2]));
    QuarterCircle(w,h,figure, ColorLine, nrLines, image, 2);
    QuarterCircle(w,h,figure, ColorLine, nrLines, image, 4);
    return image;
}

img::EasyImage Diamond(unsigned int w, unsigned int h, std::string figure, std::vector<int> ColorBG, std::vector<int>ColorLine, unsigned int nrLines){
    img::EasyImage image(w,h);
    image.clear(img::Color(ColorBG[0], ColorBG[1], ColorBG[2]));
    QuarterCircle(w/2,h/2, figure, ColorLine,nrLines, image, 1);
    QuarterCircle(w/2,h/2, figure, ColorLine,nrLines, image, 2, w/2);
    QuarterCircle(w/2,h/2, figure, ColorLine,nrLines, image, 3, w/2, h/2);
    QuarterCircle(w/2,h/2, figure, ColorLine,nrLines, image, 4, 0, h/2);
    return image;
}

img::EasyImage Lines(unsigned int w, unsigned int h, std::string figure, std::vector<int> ColorBG, std::vector<int> ColorLine, unsigned int nrLines) {
    //Lines( width, height, string figure = which kind of figure, Color of the BackGround, Color of the Lines, the amount of lines)
    img::EasyImage image(w,h);
    image.clear(img::Color(ColorBG[0], ColorBG[1], ColorBG[2]));
    if (figure == "Eye") {
        return Eye(w, h, figure, ColorBG, ColorLine, nrLines);
    }
    else if (figure == "QuarterCircle") {
        return QuarterCircle(w, h, figure, ColorLine, nrLines, image);
    }
    else if (figure == "Diamond"){
        return Diamond(w, h, figure, ColorBG, ColorLine, nrLines);
    }


}
