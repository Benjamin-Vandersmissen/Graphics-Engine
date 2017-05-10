//
// Created by Benjamin on 01/03/2017.
//

#include "WireFrameParser.hh"

WireFrameParser::WireFrameParser(const ini::Configuration &configuration, unsigned int ZBuffering) : ZBuffering(ZBuffering) {
    int size = configuration["General"]["size"].as_int_or_die();
    img::Color bgcolor = img::Color(extractColor(configuration["General"]["backgroundcolor"].as_double_tuple_or_die()));
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    int nrLights = configuration["General"]["nrLights"].as_int_or_default(0);
    bool shadowsEnabled = configuration["General"]["shadowEnabled"].as_bool_or_default(false);
    unsigned int shadowmaskSize = configuration["General"]["shadowMask"].as_int_or_default(0);
    std::vector<double> EyeCoords = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eye;
    eye = Vector3D::point(EyeCoords[0], EyeCoords[1], EyeCoords[2]);
    Figures3D figures;
    Lights3D lights;
    if (nrLights == 0){
        Light defaultAmbient;
        defaultAmbient.ambientLight = Color(1,1,1);
        lights.push_back(defaultAmbient);
    }
    for(int i = 0; i < nrLights; i++){
        Light newLight;
        std::string name = "Light" + std::to_string(i);
        newLight.ambientLight = configuration[name]["ambientLight"].as_double_tuple_or_default({0,0,0});
        newLight.infinity = configuration[name]["infinity"].as_bool_or_default(true);
        newLight.diffuseLight = configuration[name]["diffuseLight"].as_double_tuple_or_default({0,0,0});
        newLight.specularLight = configuration[name]["specularLight"].as_double_tuple_or_default({0,0,0});
        std::vector<double> vec;
        if(configuration[name]["direction"].as_double_tuple_if_exists(vec)){
            newLight.direction = Vector3D::vector(vec[0],vec[1],vec[2]);
        }else{
            newLight.direction = Vector3D::vector(0,0,0);
        }

        if (configuration[name]["location"].as_double_tuple_if_exists(vec)){
            newLight.location = Vector3D::point(vec[0],vec[1],vec[2]);
        }
        else{
            newLight.location = Vector3D::point(0,0,0);
        }
        lights.push_back(newLight);
    }

    std::map<std::string, img::EasyImage*> allTextures = {}; //asociates a path with a pointer to an image
    allTextures[""] = NULL; //default

    for(int i = 0; i < nrFigures; i++){
        std::string name = "Figure" + std::to_string(i);
        std::string type = configuration[name.c_str()]["type"].as_string_or_die();
        double rotateX = toRadial(configuration[name.c_str()]["rotateX"].as_double_or_die());
        double rotateY = toRadial(configuration[name.c_str()]["rotateY"].as_double_or_die());
        double rotateZ = toRadial(configuration[name.c_str()]["rotateZ"].as_double_or_die());
        double scale = configuration[name.c_str()]["scale"].as_double_or_die();
        std::vector<double> CenterCoords = configuration[name.c_str()]["center"].as_double_tuple_or_die();
        Vector3D CenterPoint;
        CenterPoint = Vector3D::point(CenterCoords[0], CenterCoords[1], CenterCoords[2]);
        Color ambientReflection = configuration[name.c_str()]["ambientReflection"].as_double_tuple_or_default(configuration[name.c_str()]["color"].as_double_tuple_or_default({0,0,0}));
        Color diffuseReflection = configuration[name.c_str()]["diffuseReflection"].as_double_tuple_or_default({0,0,0});
        Color specularReflection = configuration[name.c_str()]["specularReflection"].as_double_tuple_or_default({0,0,0});
        img::Color color = extractColor(configuration[name.c_str()]["color"].as_double_tuple_or_default({0,0,0}));
        double reflectionCoefficient = configuration[name.c_str()]["reflectionCoefficient"].as_double_or_default(0);
        bool rainbow = configuration[name.c_str()]["rainbow"].as_bool_or_default(false);
        std::string texturePath = configuration[name.c_str()]["texturePath"].as_string_or_default("");

        Matrix m;

        Vector3D vec;
        std::vector<Figure3D> tempfigures = {};
        if (type == "LineDrawing"){
            tempfigures.push_back(this->parseLinedrawing(configuration, name));
        }
        else if(type == "Cube"){
            tempfigures.push_back(this->parseCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if (type == "Cuboid"){
            tempfigures.push_back(this->parseCuboid(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if(type == "Tetrahedron"){
            tempfigures.push_back(this->parseTetrahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if (type == "Octahedron"){
            tempfigures.push_back(this->parseOctahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if (type == "Icosahedron"){
            tempfigures.push_back(this->parseIcosahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if (type == "Dodecahedron"){
            tempfigures.push_back(this->parseDodecahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if (type == "Cone"){
            tempfigures.push_back(this->parseCone(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if (type == "Cylinder"){
            tempfigures.push_back(this->parseCylinder(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if (type == "Sphere"){
            tempfigures.push_back(this->parseSphere(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        else if (type == "Torus"){
            tempfigures.push_back(this->parseTorus(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
//        else if (type == "Mobius"){
//            figures.push_back(this->parseMobius(configuration, name, color));
//        }
        else if (type == "3DLSystem"){
            tempfigures.push_back(this->parse3DLsystem(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
//        else if (type == "NavelTorus"){
//            figures.push_back(this->parseNavelTorus(configuration, name, color));
//        }
        else if (type == "Rail"){
            tempfigures = this->parseRail();
        }

        else if (type == "Train"){
            tempfigures = this->parseTrain(ambientReflection);
        }
        else if (type == "Station"){
            tempfigures = this->parseStation(ambientReflection);
        }

        else if (type == "Directions"){
            tempfigures = this->parseDirections();
        }
        else if (type.substr(0,7) == "Fractal"){
            tempfigures = this->parseFractal(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
        }
        else if (type.substr(0,5) == "Thick"){
            tempfigures = this->parseThick(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
        }
        else if (type == "BuckyBall"){
            tempfigures.push_back(this->parseBuckyBall(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }

        else if (type == "MengerSponge"){
            tempfigures.push_back(this->parseMengerSponge(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient));
        }
        img::Color col = {0,0,0};
        if (allTextures.find(texturePath) == allTextures.end()){
            try{
                std::ifstream textureFile(texturePath);
                img::EasyImage* newTexture = new img::EasyImage;
                textureFile >> *newTexture;
                allTextures[texturePath] = newTexture;
            }
            catch(std::exception& e){
                std::cerr << "invalid texture: " << texturePath << std::endl;
                std::cerr << e.what() << std::endl;
                allTextures[texturePath] = NULL;
            }
        }
        for(Figure3D& figure: tempfigures) {

            figure.rainbow = rainbow;
            figure.setColor(color);

            m = scaleFigure(scale);
            figure.applyTransformation(m); //apply scaling

            m = rotateFigureX(rotateX);
            figure.applyTransformation(m); //rotateFigure X-axis

            m = rotateFigureY(rotateY);
            figure.applyTransformation(m); //rotateFigure Y-axis

            m = rotateFigureZ(rotateZ);
            figure.applyTransformation(m); //rotateFigure Z-axis

            vec = Vector3D::vector(+CenterPoint.x, +CenterPoint.y,
                                   +CenterPoint.z); //translate the whole Figure if the centerpoint isn't (0,0,0);
            m = translateFigure(vec);
            figure.applyTransformation(m);
            figure.texture = allTextures[texturePath];
        }
        figures.insert(figures.end(), tempfigures.begin(), tempfigures.end());
    }
    Matrix eyeMatrix = getEyeMatrix(eye);

    applyTransformation(figures, eyeMatrix);
    Lines2D lines = doProjection(figures, ZBuffering != 0);
    img::EasyImage image;
    if (ZBuffering < 2)
        image = draw2DLines(lines, size, bgcolor, ZBuffering==1);
    else if (ZBuffering == 2){
        for(Light& light: lights){
            if (!light.infinity && shadowsEnabled){
                Matrix lightMatrix = getEyeMatrix(light.location);
                Matrix invEyeMatrix = Matrix::inv(eyeMatrix);
                applyTransformation(figures, invEyeMatrix);
                applyTransformation(figures, lightMatrix);
                Lines2D tempLines = doProjection(figures, true);
                double Xmin = tempLines.front().point1.x, Xmax = tempLines.front().point1.x, Ymin = tempLines.front().point1.y, Ymax = tempLines.front().point1.y;
                for(Line2D line: tempLines){
                    Xmin = std::min(Xmin, std::min(line.point1.x,line.point2.x));
                    Xmax = std::max(Xmax, std::max(line.point1.x,line.point2.x));
                    Ymin = std::min(Ymin, std::min(line.point1.y,line.point2.y));
                    Ymax = std::max(Ymax, std::max(line.point1.y,line.point2.y));
                }
                double Imagex = shadowmaskSize * (Xmax-Xmin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
                double Imagey = shadowmaskSize * (Ymax-Ymin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
                light.shadowmask = ZBuffer(Imagex, Imagey);
                double d = 0.95 * (Imagex/(Xmax-Xmin));
                double DCx = d * (Xmin+Xmax)/2;
                double DCy = d * (Ymin+Ymax)/2;
                double dx = (Imagex/2) - DCx;
                double dy = (Imagey/2) - DCy;
                img::EasyImage garbageImage(Imagex, Imagey);
                for (Figure3D& figure: figures){
                    triangulate(figure);
                    for(const Face& face : figure.getFaces()){
                        Vector3D pointA = figure[face[0]];
                        Vector3D pointB = figure[face[1]];
                        Vector3D pointC = figure[face[2]];
                        Lights3D emptyLights = {};
                        draw_zbuf_triangle(light.shadowmask, garbageImage, pointA, pointB, pointC, d, dx, dy,
                                           emptyLights, figure.getAmbientReflection(), figure.getDiffuseReflection(),
                                           figure.getSpecularReflection(), figure.getReflectionCoefficient(),
                                           Vector3D::point(0,0,0));
                    }
                }
                light.d = d;
                light.dx = dx;
                light.dy = dy;
                lightMatrix.inv();
                applyTransformation(figures, lightMatrix);
                applyTransformation(figures, eyeMatrix);
            }
            light.location *= eyeMatrix;
            light.direction *= eyeMatrix;
        }
        double Xmin = lines.front().point1.x, Xmax = lines.front().point1.x, Ymin = lines.front().point1.y, Ymax = lines.front().point1.y;
        for(Line2D line: lines){
            Xmin = std::min(Xmin, std::min(line.point1.x,line.point2.x));
            Xmax = std::max(Xmax, std::max(line.point1.x,line.point2.x));
            Ymin = std::min(Ymin, std::min(line.point1.y,line.point2.y));
            Ymax = std::max(Ymax, std::max(line.point1.y,line.point2.y));
        }
        double Imagex = size * (Xmax-Xmin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
        double Imagey = size * (Ymax-Ymin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
        image = img::EasyImage(roundToInt(Imagex), roundToInt(Imagey), bgcolor);
        ZBuffer buffer(roundToInt(Imagex), roundToInt(Imagey));
        double d = 0.95 * (Imagex/ (Xmax-Xmin));
        double DCx = d * (Xmin+Xmax)/2;
        double DCy = d * (Ymin+Ymax)/2;
        double dx = (Imagex/2) - DCx;
        double dy = (Imagey/2) - DCy;
        for (Figure3D& figure: figures){
            triangulate(figure);
            for(const Face& face : figure.getFaces()){
                Vector3D pointA = figure[face[0]];
                Vector3D pointB = figure[face[1]];
                Vector3D pointC = figure[face[2]];
                draw_zbuf_triangle(buffer, image, pointA, pointB, pointC, d, dx, dy,
                                   lights, figure.getAmbientReflection(), figure.getDiffuseReflection(),
                                   figure.getSpecularReflection(), figure.getReflectionCoefficient(), eye);
            }
        }
    }
    else if (ZBuffering == 3){
        double Xmin = lines.front().point1.x, Xmax = lines.front().point1.x, Ymin = lines.front().point1.y, Ymax = lines.front().point1.y;
        for(Line2D line: lines){
            Xmin = std::min(Xmin, std::min(line.point1.x,line.point2.x));
            Xmax = std::max(Xmax, std::max(line.point1.x,line.point2.x));
            Ymin = std::min(Ymin, std::min(line.point1.y,line.point2.y));
            Ymax = std::max(Ymax, std::max(line.point1.y,line.point2.y));
        }
        double Imagex = size * (Xmax-Xmin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
        double Imagey = size * (Ymax-Ymin)/(std::max((Xmax-Xmin), (Ymax-Ymin)));
        image = img::EasyImage(roundToInt(Imagex), roundToInt(Imagey), bgcolor);
        ZBuffer buffer(roundToInt(Imagex), roundToInt(Imagey));
        double d = 0.95 * (Imagex/ (Xmax-Xmin));
        double DCx = d * (Xmin+Xmax)/2;
        double DCy = d * (Ymin+Ymax)/2;
        double dx = (Imagex/2) - DCx;
        double dy = (Imagey/2) - DCy;
        Matrix invEyeMatrix = Matrix::inv(getEyeMatrix(eye));
        applyTransformation(figures, invEyeMatrix);
        unsigned int rechthoeken = 0;
        unsigned int i = 0;
        std::map<unsigned int, unsigned int> trianglesToRectangles;
        std::map<unsigned int, std::vector<Vector3D> > rectangleProperties;
        std::map<unsigned int, img::EasyImage*> rectangleTextures;
        for(Figure3D& figure: figures){
            for(Face face : figure.getFaces()){
                unsigned int rechthoek = 0;
                Vector3D A = figure[face[0]];
                Vector3D B = figure[face[1]];
                Vector3D C = figure[face[2]];
                Vector3D AB = B-A;
                Vector3D AC = C-A;
                Vector3D normal = Vector3D::normalise(Vector3D::cross(AB, AC));
                double a1 = normal.x;
                double b1 = normal.y;
                double c1 = normal.z;
                double d1 = -(a1*A.x + b1*A.y+c1*A.z);

                Vector3D linksonder;
                Vector3D linksboven;
                Vector3D rechtsonder;
                if (face.getPointIndices().size() <= 4){
                    Vector3D loodrechte = Vector3D::cross(normal, figure[face[1]]-figure[face[0]]);
                    linksonder = figure[face[0]];
                    rechtsonder = figure[face[1]];
                    linksboven = figure[face[0]]+loodrechte;
                }
                else if (face.getPointIndices().size() == 5){
                    //regelmatige vijfhoek alleen
                    Vector3D v = figure[face[face.getPointIndices().size()-1]] - figure[face[0]];
                    linksonder = figure[face[0]] - std::cos(18*M_PI/360)*v;
                    rechtsonder = figure[face[0]] + (1+std::cos(18*M_PI/360))*v;
                    Vector3D v2 = figure[face[1]] - linksonder;
                    linksboven = linksonder + (1+std::cos(36*M_PI/360))*v2;
                }
                else if (face.getPointIndices().size() == 6){
                    //regelmatige zeshoek alleen
                    Vector3D v = figure[face[face.getPointIndices().size()-1]] - figure[face[0]];
                    linksonder = figure[face[0]] - std::sin(60*M_PI/360)*v;
                    rechtsonder = figure[face[0]] + (1 + std::sin(60*M_PI/360))*v;
                    Vector3D v2 = figure[face[1]] - linksonder;
                    linksboven  = linksonder + 2*v2;
                }

                rechthoek = rechthoeken;
                //TODO: code implementeren die verschillende faces met hetzelfde vlak van een figuur dezelfde rechthoek toewijst
//                for(int j = 0; j <rectangleProperties.size();j++){
//                    Vector3D P2 = rectangleProperties[j][0];
//                    Vector3D A2 = rectangleProperties[j][1];
//                    Vector3D B2 = rectangleProperties[j][2];
//                    Vector3D normal1 = Vector3D::cross(P2+A2, P2+B2);
//                    normal1.normalise();
//                    double a2,b2,c2,d2;
//                    a2 = normal.x;
//                    b2 = normal.y;
//                    c2 = normal.z;
//                    d2 = -(a2*P2.x+b2*P2.y+c2*P2.z);
//                    if (std::abs(a2-a1) <= std::pow(10,-10) && std::abs(b1-b2) <= std::pow(10,-10) && std::abs(c1-c2) <= std::pow(10,-10) && std::abs(d1-d2) <= std::pow(10,-10)){
//                        Vector3D P1 = linksonder;
//                        Vector3D A1 = rechtsonder - linksonder;
//                        Vector3D B1 = linksboven - linksonder;
//                        Vector3D P3 = P2+A2;
//                        Vector3D P4 = P2+B2;
//                        double u1,v1,u2,v2,u3,v3;
//                        if (std::abs(A1.x*B1.y-A1.y*B1.x) >= std::pow(10,-10)){
//                            double det = (A1.x*B1.y-A1.y*B1.x);
//                            u1 = (1/det)*((P2.x-P1.x)*B1.y + (P2.y-P1.y)*(-B1.x));
//                            v1 = (1/det)*((P2.x-P1.x)*(-A1.y) + (P2.y-P1.y)*A1.x);
//                            u2 = (1/det)*((P3.x-P1.x)*B1.y + (P3.y-P1.y)*(-B1.x));
//                            v2 = (1/det)*((P3.x-P1.x)*(-A1.y) + (P3.y-P1.y)*A1.x);
//                            u3 = (1/det)*((P4.x-P1.x)*B1.y + (P4.y-P1.y)*(-B1.x));
//                            v3 = (1/det)*((P4.x-P1.x)*(-A1.y) + (P4.y-P1.y)*A1.x);
////                        std::cerr << "U1, V1: " << u << ' ' << v << std::endl;
//                        }
//                        if (std::abs(A1.y*B1.z-A1.z*B1.y) >= std::pow(10,-10)){
//                            u1 = (1/(A1.y*B1.z-A1.z*B1.y))*((P2.y-P1.y)*B1.z + (P2.z-P1.z)*(-B1.y));
//                            v1 = (1/(A1.y*B1.z-A1.z*B1.y))*((P2.y-P1.y)*(-A1.z) + (P2.z-P1.z)*A1.y);
//                            u2 = (1/(A1.y*B1.z-A1.z*B1.y))*((P3.y-P1.y)*B1.z + (P3.z-P1.z)*(-B1.y));
//                            v2 = (1/(A1.y*B1.z-A1.z*B1.y))*((P3.y-P1.y)*(-A1.z) + (P3.z-P1.z)*A1.y);
//                            u3 = (1/(A1.y*B1.z-A1.z*B1.y))*((P4.y-P1.y)*B1.z + (P4.z-P1.z)*(-B1.y));
//                            v3 = (1/(A1.y*B1.z-A1.z*B1.y))*((P4.y-P1.y)*(-A1.z) + (P4.z-P1.z)*A1.y);
////                        std::cerr << "U2, V2: " << u << ' ' << v << std::endl;
//                        }
//                        if (std::abs(A1.x*B1.z-A1.z*B1.x) >= std::pow(10,-10)){
//                            u1 = (1/(A1.x*B1.z-A1.z*B1.x))*((P2.x-P1.x)*B1.z + (P2.z-P1.z)*(-B1.x));
//                            v1 = (1/(A1.x*B1.z-A1.z*B1.x))*((P2.x-P1.x)*(-A1.z) + (P2.z-P1.z)*A1.x);
//                            u2 = (1/(A1.x*B1.z-A1.z*B1.x))*((P3.x-P1.x)*B1.z + (P3.z-P1.z)*(-B1.x));
//                            v2 = (1/(A1.x*B1.z-A1.z*B1.x))*((P3.x-P1.x)*(-A1.z) + (P3.z-P1.z)*A1.x);
//                            u3 = (1/(A1.x*B1.z-A1.z*B1.x))*((P4.x-P1.x)*B1.z + (P4.z-P1.z)*(-B1.x));
//                            v3 = (1/(A1.x*B1.z-A1.z*B1.x))*((P4.x-P1.x)*(-A1.z) + (P4.z-P1.z)*A1.x);
////                        std::cerr << "U3, V3: " << u << ' ' << v << std::endl;
//                        }
//
//                        std::cerr << rectangleProperties[j][0] << std::endl;
//                        std::cerr << rectangleProperties[j][1] << std::endl;
//                        std::cerr << rectangleProperties[j][2] << std::endl;
//                        if (u1 < 0 && v1 < 0){
//                            rectangleProperties[j][0] = P1+ u1*A1 +v1*B1;
//                        }
//                        else if (u1 > 0 && v1 < 0){
//                            rectangleProperties[j][0] = P1  + v1*B1;
//                        }
//                        else if (u1 < 0 && v1 >0){
//                            rectangleProperties[j][0] = P1 + u1*A1 ;
//                        }
//                        if (u2 > 1) {
//                            if (v2 > 0)
//                                rectangleProperties[j][1] = (P1 + u2 * A1) - rectangleProperties[j][0];
//                            else
//                                rectangleProperties[j][1] = (P1 + u2*A1 +v2*B1) - rectangleProperties[j][0];
//                        }
//                        else {
//                            if (v2 >0)
//                                rectangleProperties[j][1] = (P1 + A1) - rectangleProperties[j][0];
//                            else
//                                rectangleProperties[j][1] = (P1 + A1 + v2*B1) - rectangleProperties[j][0];
//                        }
//                        if (v3 > 1){
//                            if (u3 > 0)
//                                rectangleProperties[j][2] = (P1 + v3 * B1) - rectangleProperties[j][0];
//                            else
//                                rectangleProperties[j][2] = (P1 + u3 * A1 + v3 * B1) - rectangleProperties[j][0];
//                        }
//                        else{
//                            if (u3 > 0)
//                                rectangleProperties[j][2] = (P1 + B1) - rectangleProperties[j][0];
//                            else
//                                rectangleProperties[j][2] = (P1 + u3*A1 + B1) - rectangleProperties[j][0];
//                        }
//                        std::cerr << rectangleProperties[j][0] << ' ' << P1 << std::endl;
//                        std::cerr << rectangleProperties[j][1] << ' ' << A1 << std::endl;
//                        std::cerr << rectangleProperties[j][2] << ' ' << B1 << std::endl;
//                        rechthoek = j;
//                        break;
//
//                    }
//                }
                for(int j = 2; j < face.getPointIndices().size(); j++){ //met het trianguleren wordt een face met n hoekpunten opgesplitst in n-2 driehoeken
                    trianglesToRectangles[i] = rechthoek;
                    i++;
                }
//                std::cerr << linksonder << ' ' << rechtsonder-linksonder << ' ' << linksboven - linksonder << std::endl;


                if (rechthoek == rechthoeken) {
                    rectangleProperties[rechthoeken] = {linksonder, rechtsonder-linksonder, linksboven-linksonder};
                    rectangleTextures[rechthoeken] = figure.texture;
                    rechthoeken++;
                }
            }
        }
        applyTransformation(figures, eyeMatrix);
        int triangle = 0;
        for (Figure3D& figure: figures){
            triangulate(figure);
            for(const Face& face : figure.getFaces()){
                Vector3D pointA = figure[face[0]];
                Vector3D pointB = figure[face[1]];
                Vector3D pointC = figure[face[2]];
                unsigned int rectangle = trianglesToRectangles[triangle];
                std::vector<Vector3D > rectangleProps = rectangleProperties[rectangle];
                draw_textured_triangle(buffer, image, pointA, pointB, pointC, d, dx, dy,
                                       rectangleProps, figure.texture, eye);
                triangle++;
            }
        }
    }
    this->image = image;
}

Figure3D WireFrameParser::parseLinedrawing(const ini::Configuration &configuration, std::string &name) {
    int nrPoints = configuration[name.c_str()]["nrPoints"].as_int_or_die();
    int nrLines = configuration[name.c_str()]["nrLines"].as_int_or_die();
    std::string str;
    std::vector<Vector3D> Points;
    for(int i = 0; i < nrPoints; i++){
        str = "point" + std::to_string(i);
        std::vector<double> Coordinates = configuration[name.c_str()][str.c_str()].as_double_tuple_or_die();
        Vector3D point;
        point = Vector3D::point(Coordinates[0], Coordinates[1], Coordinates[2]);
        Points.push_back(point);
    }
    std::vector<Face> Faces;
    for(int i = 0; i < nrLines; i++){
        str = "line" + std::to_string(i);
        std::vector<int> indices = configuration[name.c_str()][str.c_str()].as_int_tuple_or_die();
        for(int j = 0; j < indices.size(); j++){
            if (indices[j] >= Points.size() || indices[j] < 0){
                std::cerr << "Index out of bounds" << std::endl;
                throw(std::invalid_argument("Index out of bounds"));
            }
        }
        Faces.push_back(Face(indices));
    }
    Figure3D figure;
    figure.setFaces(Faces);
    figure.setPoints(Points);
    return figure;
}

const img::EasyImage &WireFrameParser::getImage() const {
    return image;
}

Figure3D WireFrameParser::parseCube(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                    unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};
    return this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D WireFrameParser::parseTetrahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                           unsigned int reflectionCoefficient) {
    return this->drawTetrahedron(ambientReflection, diffuseReflection, specularReflection,reflectionCoefficient);
}

Figure3D WireFrameParser::parseOctahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                          unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};

    return this->drawOctahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D WireFrameParser::parseIcosahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                           unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};

    return this->drawIcosahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D WireFrameParser::parseDodecahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                            unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};

    return this->drawDodecahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D
WireFrameParser::parseCone(const ini::Configuration &configuration, std::string &name, Color ambientReflection, Color diffuseReflection,
                           Color specularReflection, unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};

    double height = configuration[name]["height"].as_double_or_die();
    int n = configuration[name]["n"].as_int_or_die();
    return this->drawCone(height, n, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D WireFrameParser::parseCuboid(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                                      Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};

    double height = configuration[name]["height"].as_double_or_die();
    double length = configuration[name]["length"].as_double_or_die();
    double depth = configuration[name]["depth"].as_double_or_die();
    return this->drawCuboid(length, height, depth, ambientReflection, diffuseReflection, specularReflection,
                            reflectionCoefficient, Vector3D(), Vector3D(), 0);
}

Figure3D WireFrameParser::parseCylinder(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                                        Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};

    double height = configuration[name]["height"].as_double_or_die();
    int n = 2*configuration[name]["n"].as_int_or_die();
    return this->drawCylinder(height, n, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D WireFrameParser::parseSphere(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                                      Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};

    int n = configuration[name]["n"].as_int_or_die();
    return this->drawSphere(n, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D
WireFrameParser::parseTorus(const ini::Configuration &configuration, std::string &name, Color ambientReflection, Color diffuseReflection,
                            Color specularReflection, unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};

    double r = configuration[name]["r"].as_double_or_die();
    double R = configuration[name]["R"].as_double_or_die();
    int m = configuration[name]["m"].as_int_or_die();
    int n = configuration[name]["n"].as_int_or_die();
    return this->drawTorus(r, R, m, n, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D
WireFrameParser::parse3DLsystem(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                                Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient) {
    img::Color color = ambientReflection.toRGB();

    std::string filename = configuration[name]["inputfile"].as_string_or_die();
    LParser::LSystem3D Lsystem = getLSystem3D(filename);
    Figure3D figure = drawLSystem3D(Lsystem, color);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    figure.setSpecularReflection(specularReflection);
    return figure;
}

Figure3D WireFrameParser::drawCube(Color ambientReflection, Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient,
                                   Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure;
    std::vector<Vector3D> points= {Vector3D::point(1,-1,-1), Vector3D::point(-1,1,-1), Vector3D::point(1,1,1), Vector3D::point(-1,-1,1),
                                   Vector3D::point(1,1,-1), Vector3D::point(-1,-1,-1), Vector3D::point(1,-1,1), Vector3D::point(-1,1,1)};
    figure.setPoints(points);
    std::vector<Face> faces = {Face({0,4,2,6}),Face({4,1,7,2}), Face({1,5,3,7}), Face({5,0,6,3}), Face({6,2,7,3}), Face({0,5,1,4})};
    figure.setFaces(faces);
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D WireFrameParser::drawTetrahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                          unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(1,-1,-1), Vector3D::point(-1,1,-1), Vector3D::point(1,1,1), Vector3D::point(-1,-1,1)};
    figure.setPoints(points);
    std::vector<Face> faces = {Face({0,1,2}), Face({1,3,2}), Face({0,3,1}), Face({0,2,3})};
    figure.setFaces(faces);
    
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D WireFrameParser::drawOctahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                         unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(1,0,0), Vector3D::point(0,1,0), Vector3D::point(-1,0,0), Vector3D::point(0,-1,0),
                                    Vector3D::point(0,0,-1), Vector3D::point(0,0,1)};
    figure.setPoints(points);
    std::vector<Face> faces = {Face({0,1,5}), Face({1,2,5}), Face({2,3,5}), Face({3,0,5}),
                               Face({1,0,4}), Face({2,1,4}), Face({3,2,4}), Face({0,3,4})};
    figure.setFaces(faces);
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D WireFrameParser::drawIcosahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                          unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure;
    std::vector<Vector3D> points = {Vector3D::point(0,0, std::sqrt(5)/2)};
    for(int i = 2; i < 12; i ++){
        Vector3D point;
        if(i <= 6){
            point = point.point(cos((i-2)*2*M_PI/5), sin((i-2)*2*M_PI/5), 0.5);
        }
        else{
            point = point.point(cos(M_PI/5 + (i-7)*2*M_PI/5), sin(M_PI/5 +(i-7)*2*M_PI/5), -0.5);
        }
        points.push_back(point);
    }
    points.push_back(Vector3D::point(0,0,-sqrt(5)/2));
    figure.setPoints(points);
    std::vector<Face> faces = {Face({0,1,2}), Face({0,2,3}), Face({0,3,4}), Face({0,4,5}), Face({0,5,1}),
                               Face({1,6,2}), Face({2,6,7}), Face({2,7,3}), Face({3,7,8}), Face({3,8,4}),
                               Face({4,8,9}), Face({4,9,5}), Face({5,9,10}), Face({5,10,1}), Face({1,10,6}),
                               Face({11,7,6}), Face({11,8,7}), Face({11,9,8}), Face({11,10,9}), Face({11,6,10})};
    figure.setFaces(faces);
    
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D WireFrameParser::drawDodecahedron(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                           unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure = parseIcosahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
    Figure3D figure2;
    std::vector<Vector3D> points = {};
    for(int i = 0; i < figure.getFaces().size();i++){
        points.push_back(figure.getCenter(i));
    }
    figure2.setPoints(points);
    std::vector<Face> faces = {Face({0,1,2,3,4}), Face({0,5,6,7,1}), Face({1,7,8,9,2}), Face({2,9,10,11,3}),
                               Face({3,11,12,13,4}), Face({4,13,14,5,0}), Face({19,18,17,16,15}),
                               Face({19,14,13,12,18}), Face({18,12,11,10,17}), Face({17,10,9,8,16}),
                               Face({16,8,7,6,15}), Face({15,6,5,14,19})};

    figure2.setFaces(faces);
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure2.setAmbientReflection(ambientReflection);
    figure2.setDiffuseReflection(diffuseReflection);
    figure2.setSpecularReflection(specularReflection);
    figure2.setReflectionCoefficient(reflectionCoefficient);
    return figure2;
}

Figure3D
WireFrameParser::drawCone(double height, int n, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                          unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure;
    Vector3D top = Vector3D::point(0,0,height);
    std::vector<Vector3D> points = {};
    std::vector<Face> faces = {};
    std::vector<int> circleIndices = {};
    for(int i = 0; i < n; i++){
        points.push_back(Vector3D::point(cos(2*i*M_PI/n), sin(2*i*M_PI/n) , 0));
        faces.push_back(Face({i, (i+1)%n, n}));
        circleIndices.insert(circleIndices.begin(), i);
    }
    faces.push_back(Face(circleIndices));
    points.push_back(top);
    
    figure.setFaces(faces);
    figure.setPoints(points);
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D
WireFrameParser::drawCuboid(double length, double height, double depth, Color ambientReflection, Color diffuseReflection,
                            Color specularReflection, unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation,
                            double scale) {
    Figure3D figure;
    std::vector<Vector3D> points= {Vector3D::point(length/2,-height/2,-depth/2), Vector3D::point(-length/2,height/2,-depth/2), Vector3D::point(length/2,height/2,depth/2), Vector3D::point(-length/2,-height/2,depth/2),
                                   Vector3D::point(length/2,height/2,-depth/2), Vector3D::point(-length/2,-height/2,-depth/2), Vector3D::point(length/2,-height/2,depth/2), Vector3D::point(-length/2,height/2,depth/2)};
    figure.setPoints(points);
    
    std::vector<Face> faces = {Face({0,4,2,6}),Face({4,1,7,2}), Face({1,5,3,7}), Face({5,0,6,3}), Face({6,2,7,3}), Face({0,5,1,4})};
    figure.setFaces(faces);
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D WireFrameParser::drawCylinder(double height, int n, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                       unsigned int reflectionCoefficient, bool zijvlakken, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure;
    std::vector<Vector3D> points = {};
    std::vector<Face> faces = {};
    std::vector<int> circleIndices = {};

    n*=2;
    for(int i = 0; i < n; i ++){
        if (i < n/2) {
            points.push_back(Vector3D::point(cos(4 * i * M_PI / n), sin(4 * i * M_PI / n), 0));
            faces.push_back(Face({ (i+1)%(n/2), (i+1)%(n/2) + n/2 , i + n/2, i}));
            circleIndices.push_back(i);
        }
        if (i == n/2){
            faces.push_back(Face(circleIndices));
            circleIndices = {};
        };
        if(i >= n/2){
            points.push_back(Vector3D::point(cos(4 * i * M_PI / n), sin(4 * i * M_PI / n), height));
            circleIndices.push_back(i);
        }

    }
    if (zijvlakken)
        faces.push_back(Face(circleIndices));
    
    figure.setPoints(points);
    figure.setFaces(faces);
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D WireFrameParser::drawSphere(int n, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                     unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure = drawIcosahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
    for (int j = 0; j < n; j++) { //repeat the process n times
        std::vector<Vector3D> test = {};
        std::vector<Face> facesTest = {};
        for (Face face1 : figure.getFaces()) {
            Vector3D punt1 = figure.getPoints()[face1.getPointIndices()[0]];
            Vector3D punt2 = figure.getPoints()[face1.getPointIndices()[1]];
            Vector3D punt3 = figure.getPoints()[face1.getPointIndices()[2]];
            Vector3D punt4 = (punt1 + punt2) / 2;
            Vector3D punt5 = (punt2 + punt3) / 2;
            Vector3D punt6 = (punt1 + punt3) / 2;
            int index1, index2, index3, index4, index5, index6;
            if (std::find(test.begin(), test.end(), punt1) == test.end()) {
                index1 = test.size();
                test.push_back(punt1);
            } else {
                index1 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt1));
            }
            if (std::find(test.begin(), test.end(), punt2) == test.end()) {
                index2 = test.size();
                test.push_back(punt2);
            } else {
                index2 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt2));
            }
            if (std::find(test.begin(), test.end(), punt3) == test.end()) {
                index3 = test.size();
                test.push_back(punt3);
            } else {
                index3 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt3));
            }
            if (std::find(test.begin(), test.end(), punt4) == test.end()) {
                index4 = test.size();
                test.push_back(punt4);
            } else {
                index4 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt4));
            }
            if (std::find(test.begin(), test.end(), punt5) == test.end()) {
                index5 = test.size();
                test.push_back(punt5);
            } else {
                index5 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt5));
            }
            if (std::find(test.begin(), test.end(), punt6) == test.end()) {
                index6 = test.size();
                test.push_back(punt6);
            } else {
                index6 = std::distance(test.begin(), std::find(test.begin(), test.end(), punt6));
            }
            facesTest.push_back(Face({index1, index4, index6}));
            facesTest.push_back(Face({index4, index2, index5}));
            facesTest.push_back(Face({index5, index3, index6}));
            facesTest.push_back(Face({index4, index5, index6}));
        }
        figure.setPoints(test);
        figure.setFaces(facesTest);
    }
    std::vector<Vector3D> points = figure.getPoints();
    for (Vector3D &point : points) {
        point /= point.length();
        point.normalise();
    }
    figure.setPoints(points);
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center); 
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D
WireFrameParser::drawTorus(double r, double R, int m, int n, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                           unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure;
    std::vector<Vector3D> points = {};
    std::vector<Face> faces = {};
    for(int i = 0 ; i < n; i++){
        for(int j = 0; j < m; j++){
            double u = 2*i*M_PI/n;
            double v = 2*j*M_PI/m;
            Vector3D point = Vector3D::point((R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v));
            points.push_back(point);
            faces.push_back(Face({m*i+j, m*((i+1)%n) + j, m*((i+1)%n) + (j+1)%m, m*i + (j+1)%m}));
        }
    }
    
    figure.setFaces(faces);
    figure.setPoints(points);
    Matrix matrix;
    matrix = scaleFigure(scale);
    figure.applyTransformation(matrix);
    matrix = rotateFigureX(rotation.x)*rotateFigureY(rotation.y)*rotateFigureZ(rotation.z)*translateFigure(center);
    figure.applyTransformation(matrix);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

Figure3D WireFrameParser::drawBuckyBall(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                        unsigned int reflectionCoefficient, Vector3D center, Vector3D rotation, double scale) {
    Figure3D figure = drawIcosahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
    std::vector<Vector3D> points = {};
    std::vector<Face> faces = {};
    std::vector<Face> triangles = {};
    std::vector<Vector3D> trianglePoints = {};
    for(const Face& face: figure.getFaces()){
        Vector3D Point1 = figure[face[0]]+(1.0/3.0)*(figure[face[1]] - figure[face[0]]);
        Vector3D Point2 = figure[face[0]]+(2.0/3.0)*(figure[face[1]] - figure[face[0]]);
        Vector3D Point3 = figure[face[1]]+(1.0/3.0)*(figure[face[2]] - figure[face[1]]);
        Vector3D Point4 = figure[face[1]]+(2.0/3.0)*(figure[face[2]] - figure[face[1]]);
        Vector3D Point5 = figure[face[2]]+(1.0/3.0)*(figure[face[0]] - figure[face[2]]);
        Vector3D Point6 = figure[face[2]]+(2.0/3.0)*(figure[face[0]] - figure[face[2]]);

        auto it = std::find(points.begin(), points.end(), Point1);
        int index1, index2, index3, index4, index5, index6, indexA, indexB, indexC;
        if (it == points.end()){
            index1 = points.size();
            points.push_back(Point1);
        }else{
            index1 = std::distance(points.begin(), it);
        }
        it = std::find(points.begin(), points.end(), Point2);
        if (it == points.end()){
            index2 = points.size();
            points.push_back(Point2);
        }
        else{
            index2 = std::distance(points.begin(), it);
        }
        it = std::find(points.begin(), points.end(), Point3);
        if (it == points.end()){
            index3 = points.size();
            points.push_back(Point3);
        }else{
            index3 = std::distance(points.begin(), it);
        }
        it = std::find(points.begin(), points.end(), Point4);
        if (it == points.end()){
            index4 = points.size();
            points.push_back(Point4);
        }else{
            index4 = std::distance(points.begin(), it);
        }
        it = std::find(points.begin(), points.end(), Point5);
        if (it == points.end()){
            index5 = points.size();
            points.push_back(Point5);
        }
        else{
            index5 = std::distance(points.begin(), it);
        }
        it = std::find(points.begin(), points.end(), Point6);
        if (it == points.end()){
            index6 = points.size();
            points.push_back(Point6);
        }else{
            index6 = std::distance(points.begin(), it);
        }
        it = std::find(trianglePoints.begin(), trianglePoints.end(), figure[face[0]]);
        if (it == trianglePoints.end()){
            indexA = trianglePoints.size();
            trianglePoints.push_back(figure[face[0]]);
        }else{
            indexA = std::distance(trianglePoints.begin(), it);
        }
        it = std::find(trianglePoints.begin(), trianglePoints.end(), figure[face[1]]);
        if (it == trianglePoints.end()){
            indexB = trianglePoints.size();
            trianglePoints.push_back(figure[face[1]]);
        }
        else{
            indexB = std::distance(trianglePoints.begin(), it);
        }
        it = std::find(trianglePoints.begin(), trianglePoints.end(), figure[face[2]]);
        if (it == trianglePoints.end()){
            indexC = trianglePoints.size();
            trianglePoints.push_back(figure[face[2]]);
        }
        else{
            indexC = std::distance(trianglePoints.begin(), it);
        }
        faces.push_back(Face({index1, index2, index3, index4, index5, index6})); //Zeshoek
        triangles.push_back(Face({indexA, index1, index6})); //Driehoek A
        triangles.push_back(Face({indexB, index3, index2})); //Driehoek B
        triangles.push_back(Face({indexC, index5, index4})); //Driehoek C
    }
    std::vector<Face> alGebruikteTriangles = {};
    for(Face& triangle : triangles){
        std::vector<Face> vijfhoekTriangles = {};
        for(Face& triangle2 : triangles){
            if (areAlmostEqual(trianglePoints[triangle[0]], trianglePoints[triangle2[0]])){
                vijfhoekTriangles.push_back(triangle2);
                alGebruikteTriangles.push_back(triangle2);
            }
        }
            rearrangeTriangles(vijfhoekTriangles, points);
        std::vector<int> indices = {};
        for(Face& tempTriangle: vijfhoekTriangles){
            int index = tempTriangle[1];
            indices.push_back(index);
        }
        faces.push_back(Face(indices));
    }
    figure.setFaces(faces);
    figure.setPoints(points);
    figure.setAmbientReflection(ambientReflection);
    figure.setDiffuseReflection(diffuseReflection);
    figure.setSpecularReflection(specularReflection);
    figure.setReflectionCoefficient(reflectionCoefficient);
    return figure;
}

std::vector<Figure3D> WireFrameParser::parseRail() {
    Figures3D figures = {};
    Color colGray = Color(0.5,0.5,0.5);
    Color colBrown = Color(0.5,0.22,0);
    figures.push_back(this->drawCylinder(20, 50, colGray, colGray, colGray, 20, true, Vector3D::point(-8, 0, 0)));
    figures.push_back(this->drawCylinder(20, 50, colGray, colGray, colGray, 20, true, Vector3D::point(8, 0, 0)));
    figures.push_back(
            this->drawCuboid(0.5, 16, 2, colBrown, colBrown, colBrown, 5, Vector3D::point(0, 0, 2),
                             Vector3D::vector(0, 0, toRadial(90))));
    figures.push_back(
            this->drawCuboid(0.5, 16, 2, colBrown, colBrown, colBrown, 5, Vector3D::point(0, 0, 7),
                             Vector3D::vector(0, 0, toRadial(90))));
    figures.push_back(
            this->drawCuboid(0.5, 16, 2, colBrown, colBrown, colBrown, 5, Vector3D::point(0, 0, 12),
                             Vector3D::vector(0, 0, toRadial(90))));
    figures.push_back(
            this->drawCuboid(0.5, 16, 2, colBrown, colBrown, colBrown, 5, Vector3D::point(0, 0, 17),
                             Vector3D::vector(0, 0, toRadial(90))));
    return figures;
}

std::vector<Figure3D> WireFrameParser::parseTrain(Color baseColor) {
    Figures3D figures = {};
    Color colDkGray = Color(0.42,0.42,0.42);
    Color colBlack = Color(0.05,0.05,0.05); //This ain't no VantaBlack
    figures.push_back(this->drawCylinder(0.8, 50, colDkGray, colDkGray, colDkGray, 20, true, Vector3D::point(-9, 2.5, 9),
                                         Vector3D::vector(0, M_PI / 2, 0), 2.5));
    figures.push_back(this->drawCylinder(0.8, 50, colDkGray, colDkGray, colDkGray, 20, true, Vector3D::point(7, 2.5, 9),
                                         Vector3D::vector(0, M_PI / 2, 0), 2.5));
    figures.push_back(this->drawCylinder(0.8, 50, colDkGray, colDkGray, colDkGray, 20, true, Vector3D::point(-9, 2.5, 3),
                                         Vector3D::vector(0, M_PI / 2, 0), 2.5));
    figures.push_back(this->drawCylinder(0.8, 50, colDkGray, colDkGray, colDkGray, 20, true, Vector3D::point(7, 2.5, 3),
                                         Vector3D::vector(0, M_PI / 2, 0), 2.5));
    figures.push_back(this->drawCylinder(0.8, 50, colDkGray, colDkGray, colDkGray, 20, true, Vector3D::point(-9, 2.5, 15),
                                         Vector3D::vector(0, M_PI / 2, 0), 2.5));
    figures.push_back(this->drawCylinder(0.8, 50, colDkGray, colDkGray, colDkGray, 20, true, Vector3D::point(7, 2.5, 15),
                                         Vector3D::vector(0, M_PI / 2, 0), 2.5));
    figures.push_back(this->drawCylinder(0.8, 50, colDkGray, colDkGray, colDkGray, 20, true, Vector3D::point(-9, 2.5, 21),
                                         Vector3D::vector(0, M_PI / 2, 0), 2.5));
    figures.push_back(this->drawCylinder(0.8, 50, colDkGray, colDkGray, colDkGray, 20, true, Vector3D::point(7, 2.5, 21),
                                         Vector3D::vector(0, M_PI / 2, 0), 2.5));
    figures.push_back(this->drawCuboid(16, 10, 25, baseColor, baseColor, baseColor, 15, Vector3D::point(0, 7.5, 12.5)));
    figures.push_back(this->drawCuboid(14, 4, 23, colBlack, colBlack, colBlack, 15, Vector3D::point(0, 14.5, 12.5)));
    return figures;
}

std::vector<Figure3D> WireFrameParser::parseStation(Color baseColor) {
    Figures3D figures = {};
    Color colDKGray = Color(0.42,0.42,0.42);
    figures.push_back(this->drawCuboid(25, 20, 12, baseColor, baseColor, baseColor, 5, Vector3D::point(12.5,10,4)));
    figures.push_back(this->drawCuboid(25 ,1,20, colDKGray, colDKGray, colDKGray, 20, Vector3D::point(12.5,20.5,0)));
    figures.push_back(this->drawCylinder(20,20, colDKGray, colDKGray, colDKGray, 20, true, Vector3D::point(1,0,-9), Vector3D::vector(toRadial(-90),0, 0)));
    figures.push_back(this->drawCylinder(20,20, colDKGray, colDKGray, colDKGray, 20, true, Vector3D::point(24,0,-9), Vector3D::vector(toRadial(-90),0, 0)));
    return figures;
}

std::vector<Figure3D> WireFrameParser::parseDirections() {
    Figures3D figures = {};
    std::vector<Face> faces = {Face({0,1})};
    std::vector<Vector3D> points = {Vector3D::point(0,0,0), Vector3D::point(100,0,0)};
    img::Color color = img::Color(255,0,0);
    figures.push_back(Figure3D(faces, points,color));
    color = img::Color(0,255,0);
    points[1] = Vector3D::point(0,100,0);
    figures.push_back(Figure3D(faces,points,color));
    color = img::Color(0,0,255);
    points[1] = Vector3D::point(0,0,100);
    figures.push_back(Figure3D(faces,points,color));
    return figures;
}

Figures3D WireFrameParser::parseFractal(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                                        Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient) {
    std::vector<std::string> supportedtypes = {"Cube", "Dodecahedron", "Icosahedron", "Octahedron", "Tetrahedron", "BuckyBall"};
    std::string type = configuration[name]["type"].as_string_or_die();
    type = type.substr(type.find("Fractal") + 7, type.length());
    Figure3D figure;
    img::Color color = {0,0,0};

    switch(std::distance(supportedtypes.begin(),std::find(supportedtypes.begin(), supportedtypes.end(), type))){
        case 0 :
            figure = this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 1 :
            figure = this->drawDodecahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 2 :
            figure = this->drawIcosahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 3 :
            figure = this->drawOctahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 4 :
            figure = this->drawTetrahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 5 :
            figure = this->drawBuckyBall(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        default:
            std::cerr << "Unknown type: " << configuration[name][type].as_string_or_die() << std::endl;
            return {};
    }
    const int iterations = configuration[name]["nrIterations"].as_int_or_die();
    const double fractalScale = configuration[name]["fractalScale"].as_double_or_die();
    return this->generateFractal(figure, iterations, fractalScale);
}

Figures3D WireFrameParser::generateFractal(Figure3D &figure, const int iterations, const double scale) {
    Figures3D figures = {figure};
    for(int i = 0; i < iterations; i ++){
        Figures3D newFigures = {};
        for(Figure3D& figure1: figures){
            for(int j = 0; j < figure1.getPoints().size(); j ++){
                Figure3D tempfigure = figure1;
                Matrix matrix = scaleFigure(1/scale);
                tempfigure.applyTransformation(matrix);
                Vector3D vector = figure1.getPoints()[j]-tempfigure.getPoints()[j];
                matrix = translateFigure(vector);
                tempfigure.applyTransformation(matrix);
                newFigures.push_back(tempfigure);
            }
        }
        figures = newFigures;
    }
    return figures;
}

Figure3D WireFrameParser::parseBuckyBall(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                         unsigned int reflectionCoefficient) {
    img::Color color = {0,0,0};
    return this->drawBuckyBall(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
}

Figure3D WireFrameParser::drawMengerSponge(Color ambientReflection, Color diffuseReflection, Color specularReflection,
                                           unsigned int reflectionCoefficient, const int iterations, Vector3D center, Vector3D rotation,
                                           double scale) {
    Figure3D sponge = this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
    Figures3D figures = {sponge};
    for(int i = 0; i < iterations; i++) {
        Figures3D temp = {};
        for(Figure3D figure: figures){
            Vector3D tempCenter = figure.getCenter();
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(2 * pow(1.0 / 3.0, i + 1) + tempCenter.x, tempCenter.y,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(-2 * pow(1.0 / 3.0, i + 1) + tempCenter.x, tempCenter.y,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(2 * pow(1.0 / 3.0, i + 1) + tempCenter.x, tempCenter.y,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(-2 * pow(1.0 / 3.0, i + 1) + tempCenter.x, tempCenter.y,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));

            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(-2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(-2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));

            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(-2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(-2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));

            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(tempCenter.x, 2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(tempCenter.x, -2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(tempCenter.x, 2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(tempCenter.x, -2 * pow(1.0 / 3.0, i + 1) + tempCenter.y,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));

            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.y, tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(-2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          -2 * pow(1.0 / 3.0, i + 1) + tempCenter.y, tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.y, tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
            temp.push_back(this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient,
                                          Vector3D::point(-2 * pow(1.0 / 3.0, i + 1) + tempCenter.x,
                                                          2 * pow(1.0 / 3.0, i + 1) + tempCenter.y, tempCenter.z),
                                          Vector3D::vector(0, 0, 0), pow(1.0 / 3, i + 1)));
        }
        figures = temp;
    }
    sponge = mergeFigures(figures);
    return sponge;
}

Figure3D
WireFrameParser::parseMengerSponge(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                                   Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient) {
    int iterations = configuration[name]["nrIterations"].as_int_or_die();

    return this->drawMengerSponge(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient, iterations);
}

Figure3D WireFrameParser::mergeFigures(Figures3D &figures) {
    std::vector<Vector3D> points = {};
    std::vector<Face> faces = {};
    for(Figure3D& figure: figures){
        for(const Vector3D& point : figure.getPoints()){
            if(std::find(points.begin(), points.end(),point ) == points.end()){
                points.push_back(point);
            }
        }
        for(const Face& face : figure.getFaces()){
            std::vector<int> newFace = {};
            for(const int i: face.getPointIndices()){
                newFace.push_back(std::distance(points.begin(),std::find(points.begin(), points.end(), figure[i])));
            }
            faces.push_back(Face(newFace));
        }
    }
    // alleen figuren met dezelfde kleur mogen gemerged worden..
    Figure3D figure;
    figure.setAmbientReflection(figures[0].getAmbientReflection());
    figure.setDiffuseReflection(figures[0].getDiffuseReflection());
    figure.setSpecularReflection(figures[0].getSpecularReflection());
    figure.setReflectionCoefficient(figures[0].getReflectionCoefficient());
    figure.setPoints(points);
    figure.setFaces(faces);
    return figure;
}

void WireFrameParser::rearrangeTriangles(std::vector<Face> &triangles, std::vector<Vector3D> points) {
    std::vector<Face> newTriangles = {triangles[0]};
    Face Triangle = triangles[0];
    for(int i = 1; i < triangles.size(); i++){
        int j = 0;
        for(Face& temp : triangles){
            j++;
            if(areAlmostEqual(points[temp[1]], points[Triangle[2]]) && std::find(newTriangles.begin(), newTriangles.end(), temp) == newTriangles.end()){
                Triangle = temp;
                newTriangles.push_back(temp);
            }
        }
    }
    triangles = newTriangles;
}

Figures3D
WireFrameParser::parseThick(const ini::Configuration &configuration, std::string &name, Color ambientReflection,
                            Color diffuseReflection, Color specularReflection, unsigned int reflectionCoefficient) {
    std::vector<std::string> supportedtypes = {"Cube", "Dodecahedron", "Icosahedron", "Octahedron", "Tetrahedron", "BuckyBall", "3DLSystem", "LineDrawing"};
    std::string type = configuration[name]["type"].as_string_or_die();
    type = type.substr(type.find("Thick") + 5, type.length());
    Figure3D figure;

    switch(std::distance(supportedtypes.begin(),std::find(supportedtypes.begin(), supportedtypes.end(), type))){
        case 0 :
            figure = this->drawCube(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 1 :
            figure = this->drawDodecahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 2 :
            figure = this->drawIcosahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 3 :
            figure = this->drawOctahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 4 :
            figure = this->drawTetrahedron(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 5 :
            figure = this->drawBuckyBall(ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        case 6:
            figure = this->parse3DLsystem(configuration, name, ambientReflection, diffuseReflection, specularReflection, reflectionCoefficient);
            break;
        default:
            std::cerr << "Unknown type: " << configuration[name]["type"].as_string_or_die() << std::endl;
            return {};
    }
    double radius = configuration[name]["radius"].as_double_or_die();
    int n = configuration[name]["n"].as_int_or_die();
    int m = configuration[name]["m"].as_int_or_die();
    return makeThicc(figure, radius, n, m);
}

Figures3D WireFrameParser::makeThicc(Figure3D &figure, const double radius, const int n, const int m) {
    Figures3D figures;
    for(Vector3D point : figure.getPoints()){
        figures.push_back(this->drawSphere(m, figure.getAmbientReflection(), figure.getDiffuseReflection(), figure.getSpecularReflection(), figure.getReflectionCoefficient(), point , Vector3D::vector(0,0,0), radius));
    }
    for(const Face& face : figure.getFaces()){
        for(int i = 0 ; i < face.getPointIndices().size(); i++){
            Vector3D p1 = figure[face[i]];
            Vector3D p2 = figure[face[(i+1)%face.getPointIndices().size()]];
            Vector3D line = p2-p1;
            double height = line.length() / radius;
            Vector3D Pr = line;
            double r = line.length();
            double theta = std::atan2(Pr.y,Pr.x);
            double phi = std::acos(Pr.z/r);
            figures.push_back(this->drawCylinder(height, n, figure.getAmbientReflection(), figure.getDiffuseReflection(), figure.getSpecularReflection(), figure.getReflectionCoefficient(),false, p1, Vector3D::vector(0,phi,theta), radius));
        }
    }
    return figures;
}









//
//Figure3D WireFrameParser::parseMobius(const ini::Configuration &configuration, std::string &name, img::Color &color) {
//    Figure3D figure;
//    double r = configuration[name]["r"].as_double_or_die();
//    double R = configuration[name]["R"].as_double_or_die();
//    int m = configuration[name]["m"].as_int_or_die();
//    int n = configuration[name]["n"].as_int_or_die();
//    std::vector<Vector3D> points;
//    std::vector<Face> faces;
//    for(int i = 0 ; i < n; i++){
//        for(int j = 0; j < m; j++){
//            double u = 2*i*M_PI/n;
//            double v = 2*j*M_PI/m;
//            Vector3D point = Vector3D::point((1+v/2*cos(u/2))*cos(u), (1+v/2*cos(u/2))*sin(u), v/2*sin(u/2));
//            points.push_back(point);
//            faces.push_back(Face({m*i+j, m*((i+1)%n) + j, m*((i+1)%n) + (j+1)%m, m*i + (j+1)%m}));
//        }
//    }
//    
//    figure.setFaces(faces);
//    figure.setPoints(points);
//    return figure;
//}
//
//Figure3D WireFrameParser::parseNavelTorus(const ini::Configuration &configuration, std::string &name, img::Color &color) {
//    Figure3D figure;
//    double r = configuration[name]["r"].as_double_or_die();
//    double R = configuration[name]["R"].as_double_or_die();
//    int m = configuration[name]["m"].as_int_or_die();
//    int n = configuration[name]["n"].as_int_or_die();
//    std::vector<Vector3D> points;
//    std::vector<Face> faces;
//    for(int i = 0 ; i < n; i++){
//        for(int j = 0; j < m; j++){
//            double u = 2*i*M_PI/n;
//            double v = 2*j*M_PI/m;
//            Vector3D point = Vector3D::point(sin(u)*(7+cos(u/3-2*v) + 2* cos(u/3+v)), cos(u)*(7+cos(u/3-2*v) + 2* cos(u/3+v)), sin(u/3-2*v)+2*sin(u/3+v));
//            points.push_back(point);
//            faces.push_back(Face({m*i+j, m*((i+1)%n) + j, m*((i+1)%n) + (j+1)%m, m*i + (j+1)%m}));
//        }
//    }
//    
//    figure.setFaces(faces);
//    figure.setPoints(points);
//    return figure;
//}