// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <random>

//include header for io
#include <iostream>
#include <fstream>
#include <iomanip>

#include <Eigen/Dense>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>

#include <common/shader.hpp>
#include <common/texture.hpp>
//#include <common/controls.hpp>

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

typedef Eigen::Array<double,1,Eigen::Dynamic,Eigen::RowMajor> ArrayXdRowMajor;
typedef Eigen::Array<double,Eigen::Dynamic,3,Eigen::RowMajor> ArrayX3dRowMajor;
typedef Eigen::Array<double,Eigen::Dynamic,2,Eigen::RowMajor> ArrayX2dRowMajor;
typedef Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ArrayXXdRowMajor;

GLFWwindow* window;

//define object shape
const int a = 3;
const int b = 50;
const int c = 3;

const int num_masses = (a+1)*(b+1)*(c+1);
const int num_cubes = a*b*c;
const int num_springs = 28*num_cubes;

//define some parameters of cubes
const double PI = 3.1415926;
const double Mass = 0.1;  //kg
const double SpringConstraint = 1000;  //N/m
const double Length = 0.1;
const double w = 4*PI;
const double CoefaRange = 0.08;
const double CoefbRange = 2*PI;

//friction coefficients for glass on glass
const double Muk = 0.8;
const double Mus = 1;

double InitialHeight = 0;
const double InitialVelocityY = 0.0;
const double TimeStep = 0.0005;  //s
const double Gravity[3] = {0,0,-9.81};

const double kc = 10000;

const double DampingCoef = 0.999;

inline double sq(double para)
{
    return para*para;
}

inline double cb(double para)
{
    return para*para*para;
}


Eigen::ArrayXi initCubeVertex(){
    Eigen::ArrayXi cubeVertex(8);
    cubeVertex <<
               0,
            1,
            a+1,
            a+2,
            (a+1)*(b+1),
            (a+1)*(b+1)+1,
            (a+1)*(b+1)+a+1,
            (a+1)*(b+1)+a+2;
    return cubeVertex;
}


Eigen::ArrayXi initTriangleVertex(){
    Eigen::ArrayXi triangleVertex(36);
    triangleVertex <<
                   0,1,2,3,1,2,
            4,5,6,7,5,6,
            1,3,7,1,5,7,
            0,2,6,0,4,6,
            7,3,2,7,6,2,
            5,1,0,5,4,0;
    return triangleVertex;
}


ArrayX2dRowMajor initBaseSpringtoMass() {
    ArrayX2dRowMajor baseSpringtoMass(28,2);
    baseSpringtoMass <<
                     //12 edges
                     0, 1,
            1, 3,
            3, 2,
            2, 0,
            4, 5,
            5, 7,
            7, 6,
            6, 4,
            0, 4,
            1, 5,
            3, 7,
            2, 6,
            //12 short diagonals
            0, 5,
            1, 7,
            3, 6,
            2, 4,
            1, 4,
            3, 5,
            2, 7,
            0, 6,
            0, 3,
            1, 2,
            4, 7,
            5, 6,
            //4 long diagonals
            0, 7,
            1, 6,
            3, 4,
            2, 5;
    return baseSpringtoMass;
}


ArrayX3dRowMajor initBaseMassPosition(){
    ArrayX3dRowMajor baseMassPosition(8,3);
    int mass_index = 0;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<2;k++){
                baseMassPosition.row(mass_index) << k,j,i;
                mass_index++;
            }
        }
    }
    baseMassPosition *= Length;
    //std::cout<<baseMassPosition<<std::endl;
    return baseMassPosition;
}


ArrayX3dRowMajor initMassPosition(){
    ArrayX3dRowMajor massPosition(num_masses,3);
    int mass_index = 0;
    for(int i=0;i<c+1;i++){
        for(int j=0;j<b+1;j++){
            for(int k=0;k<a+1;k++){
                massPosition.row(mass_index) << k,j,i;
                mass_index++;
            }
        }
    }
    massPosition *= Length;
    for (int i=0;i<num_masses;i++){
        massPosition(i,2) += InitialHeight;
    }
    //std::cout<<massPosition<<std::endl;
    return massPosition;
}


ArrayX3dRowMajor initMassVelocity(){
    ArrayX3dRowMajor massVelocity = ArrayX3dRowMajor::Zero(num_masses,3);
    for (int i=0;i<num_masses;i++){
        massVelocity(i,1) += InitialVelocityY;
    }
    return massVelocity;
}


ArrayX2dRowMajor initSpringtoMass(
        ArrayX2dRowMajor& baseSpringtoMass,
        ArrayX3dRowMajor& baseMassPosition){
    ArrayX2dRowMajor springtoMass(28*num_cubes,2);
    for(int i=0;i<28*num_cubes;i++) {
        int cube_index = i / 28;
        //std::cout<<cube_index<<std::endl;
        int spring_position = i % 28;
        int deltax = (cube_index % (a * b)) % a;
        int deltay = (cube_index % (a * b) / a);
        int deltaz = (cube_index / (a * b));
        for (int j = 0; j < 2; j++) {
            int mass_position = baseSpringtoMass(spring_position,j);
            //std::cout<<mass_position<<std::endl;
            int x0 = int(baseMassPosition(mass_position,0)/Length);
            int y0 = int(baseMassPosition(mass_position,1)/Length);
            int z0 = int(baseMassPosition(mass_position,2)/Length);
            springtoMass(i,j) = (z0 + deltaz) * (a + 1) * (b + 1) + (y0 + deltay) * (a + 1) + x0 + deltax;
            //std::cout<<"springtoMass(i,j)="<<springtoMass(i,j)<<std::endl;
        }
    }
    for(int i=0;i<28*a*b*c;i+=28){
        //std::cout<<springtoMass.row(i)<<std::endl;
    }
    return springtoMass;
}


ArrayXXdRowMajor initSpringCoefa(int num_cubes, int population_size){

    ArrayXXdRowMajor springCoefa(num_cubes,population_size);
    springCoefa = CoefaRange*ArrayXXdRowMajor::Random(num_cubes,population_size);
    return springCoefa;
}


ArrayXXdRowMajor initSpringCoefb(int num_cubes, int population_size){
    ArrayXXdRowMajor springCoefb(num_cubes,population_size);
    springCoefb = CoefbRange*ArrayXXdRowMajor::Random(num_cubes,population_size);
    //std::cout<<springCoefb<<std::endl;
    return springCoefb;
}


//rest length of springs in cube 0 at time 0
Eigen::ArrayXd initL0(){
    Eigen::ArrayXd l0(num_springs);
    for (int i=0;i<num_springs;i++){
        if ((i%28)<12)
            l0[i] = Length;
        else if ((i%28)<24)
            l0[i] = Length*sqrt(2);
        else
            l0[i] = Length*sqrt(3);
    }
    //std::cout<<l0<<std::endl;
    return l0;
}


inline void applyGravity(ArrayX3dRowMajor& massForces){
    ArrayX3dRowMajor all_gravity(num_masses,3);
    for (int i=0;i<num_masses;i++){
        all_gravity.row(i)<<Gravity[0],Gravity[1],Mass*Gravity[2];
    }
    massForces += all_gravity;
    //std::cout<<massForces<<std::endl;
}


inline void applySpringForces(
        Eigen::ArrayXd& l0t,
        ArrayX3dRowMajor& massPosition,
        ArrayX3dRowMajor& massForces,
        ArrayX2dRowMajor& springtoMass){

    //std::cout<<"l0t"<<l0t<<std::endl;
    for (int i=0;i<num_springs;i++){
        Eigen::Array<double,1,3,Eigen::RowMajor> lxyz;
        lxyz = (massPosition.row(springtoMass(i,0)) - massPosition.row(springtoMass(i,1)));
        double lt = lxyz.matrix().norm();

        auto force = SpringConstraint*(lt-l0t(i))*lxyz/lt;
        //std::cout<<"lt spring"<<i<<" = "<<lt<<std::endl;
        //std::cout<<"force"<<force<<std::endl;

        massForces.row(springtoMass(i,0)) = massForces.row(springtoMass(i,0)) - force;
        massForces.row(springtoMass(i,1)) = massForces.row(springtoMass(i,1)) + force;
    }
}


inline void applyFrictionandGroundForces(
        ArrayX3dRowMajor& massPosition,
        ArrayX3dRowMajor& massForces,
        ArrayX3dRowMajor& massVelocity
        ){

    for(int i=0;i<num_masses;i++){
        double height = massPosition(i,2);
        double ground_force = kc * sq(height);
        if (height<0){
            massForces(i,2) += ground_force;
            double v_horizontal = std::sqrt(sq(massVelocity(i,0))+sq(massVelocity(i,1)));


            double force_h = std::sqrt(sq(massForces(i,0))+sq(massForces(i,1)));
            if(v_horizontal>1e-15){
                // dynamic friction
                double dynamic_friction = Muk * ground_force;
                massForces(i,0) -= dynamic_friction*massVelocity(i,0)/v_horizontal;
                massForces(i,1) -= dynamic_friction*massVelocity(i,1)/v_horizontal;
            }
            else{
                double max_friction = Mus * ground_force;
                if (max_friction>=force_h) {
                    massForces(i,0) = 0;
                    massForces(i,1) = 0;
                }
                else{
                    double dynamic_friction = Muk * ground_force;
                    massForces(i,0) -= dynamic_friction*massForces(i,0)/force_h;
                    massForces(i,1) -= dynamic_friction*massForces(i,1)/force_h;
                }

            }

        }
    }
}


Eigen::ArrayXd Setl0t(
        double& timeStamp,
        ArrayXdRowMajor& springCoefa,
        ArrayXdRowMajor& springCoefb,
        Eigen::ArrayXd& l0){
    Eigen::ArrayXd l0t(num_springs);
    //std::cout<<"l0t.size()="<<l0t.size()<<std::endl;
    //std::cout<<"l0.size()="<<l0.size()<<std::endl;
    //std::cout<<"springCoefa.size()="<<springCoefa.size()<<std::endl;
    //std::cout<<"springCoefb.size()="<<springCoefb.size()<<std::endl;
    //std::cout<<"l0.size()="<<l0.size()<<std::endl;
    //std::cout<<"l0\n"<<l0<<std::endl;
    for (int i=0;i<num_springs;i++){
        l0t(i) = l0(i)+springCoefa(i/28)*sin(w*timeStamp+springCoefb(i/28));
        //std::cout<<"l0t"<<i<<"="<<l0t(i)<<std::endl;
    }
    return l0t;
}

void calculateEnergy(
        ArrayX3dRowMajor& massVelocity,
        ArrayX3dRowMajor& massPosition,
        ArrayX2dRowMajor& springtoMass,
        Eigen::ArrayXd& l0t){
    double kinetic_energy = 0;
    for (int i=0;i<num_masses;i++){
        kinetic_energy += Mass*(sq(massVelocity(i,0))+sq(massVelocity(i,1))+sq(massVelocity(i,2)))/2;
    }
    //std::cout<<"velocity = "<<gravitational_energy<<std::endl;

    double gravitational_energy = 0;
    for (int i=0;i<num_masses;i++){
        gravitational_energy -= Mass*Gravity[2]*massPosition(i,2);
        //std::cout<<"massPosition(i,2) ="<<massPosition(i,2)<<std::endl;
        //std::cout<<"gravitational i ="<<gravitational_energy<<std::endl;
    }
//    std::cout<<"height"<<massPosition(0,2)<<std::endl;

    double spring_energy = 0;
    for (int i=0;i<num_springs;i++){
        Eigen::Array<double,1,3,Eigen::RowMajor> lxyz;
        lxyz << (massPosition.row(springtoMass(i,0)) - massPosition.row(springtoMass(i,1)));
        double lt = std::sqrt(lxyz(0)*lxyz(0)+lxyz(1)*lxyz(1)+lxyz(2)*lxyz(2));
        spring_energy += SpringConstraint * sq(lt-l0t(i))/2;
    }
    //std::cout<<"spring"<<spring_energy<<std::endl;
    //std::cout<<"gravitational+kinetic+spring = "<<kinetic_energy+gravitational_energy+spring_energy<<std::endl;

    double ground_energy = 0;
    for (int i=0;i<num_masses;i++) {
        if (massPosition(i,2)<0){
            ground_energy -= kc*cb(massPosition(i,2))/3;
        }
    }
    double total_energy = kinetic_energy + gravitational_energy + spring_energy + ground_energy;
    //std::cout<<"total = "<<total_energy<<std::endl;
}


glm::vec3 cameraPos   = glm::vec3(4.0f, 12.0f,  2.0f);
glm::vec3 cameraFront = glm::vec3(-5.0f, -3.0f, -2.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 0.0f,  1.0f);


void processInput(GLFWwindow *window)
{
    float cameraSpeed = 0.05f; // adjust accordingly
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos += cameraSpeed * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos -= cameraSpeed * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;

}

void framebufferSizeCallback(GLFWwindow *window, int w, int h){
    glViewport(0,0,w,h);
}

int render(ArrayX3dRowMajor& position_history,
           Eigen::ArrayXi& triangleVertex,
           Eigen::ArrayXi& cubeVertex,
           ArrayX2dRowMajor& springtoMass){
    // Initialise GLFW
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        getchar();
        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Open a window and create its OpenGL context
    window = glfwCreateWindow(1440, 1080, "3 * 3 * 50 Cubes", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        getchar();
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        getchar();
        glfwTerminate();
        return -1;
    }

    glfwSetFramebufferSizeCallback(window,framebufferSizeCallback);
    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

    // Dark blue background
    //glClearColor(0, 0, 0.4, 0.0f);
    //or green color
    glClearColor(0.2, 0.8, 0.2, 0.0f);


    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);




    // Create and compile our GLSL program from the shaders
    GLuint programID = LoadShaders( "TransformVertexShader.vertexshader", "ColorFragmentShader.fragmentshader" );

    // Get a handle for our "MVP" uniform
    GLuint MatrixID = glGetUniformLocation(programID, "MVP");

    // Projection matrix : 45� Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
    glm::mat4 Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
    // Camera matrix




    glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);

//    glm::mat4 view       = glm::lookAt(
//            //glm::vec3(4,3,-3), // Camera is at (4,3,-3), in World Space
//            glm::vec3(5,5,2),
//            glm::vec3(0,0,0), // and looks at the origin
//            glm::vec3(0,0,1)  // Head is up (set to 0,-1,0 to look upside-down)
//    );
    // Model matrix : an identity matrix (model will be at the origin)
    glm::mat4 Model      = glm::mat4(1.0f);
    // Our ModelViewProjection : multiplication of our 3 matrices
    glm::mat4 MVP        = Projection * view * Model; // Remember, matrix multiplication is the other way around

    int k = 1;
    // Our vertices. Tree consecutive floats give a 3D vertex; Three consecutive vertices give a triangle.
    // A cube has 6 faces with 2 triangles each, so this makes 6*2=12 triangles, and 12*3 vertices
    int num_vertices_per_cube = 36;
    GLfloat g_vertex_buffer_data[3*num_vertices_per_cube*num_cubes];

    Eigen::Array<float,Eigen::Dynamic,3,Eigen::RowMajor> massPositionFloat(num_masses,3);






//    int num_triangles_per_cube = 12;
//    int num_triangles = num_triangles_per_cube*num_cubes;
//
//    Eigen::Array<unsigned int,Eigen::Dynamic,3,Eigen::RowMajor> indices(num_triangles,3);
//
//
//    for(int i=0;i<num_cubes;i++){
//        for (int j = 0; j < num_triangles_per_cube; ++j) {
//
//        }
//    }
//
//
//    //draw 36 vertices
//    for (int i=0;i<num_triangles;i++){
//        indices.row(i) << triangleVertex(
//        cubeVertex.row(triangleVertex(i))+k)
//
//        g_vertex_buffer_data[3*i+0] = massPositionFloat(cubeVertex(triangleVertex(i))+k,0);
//        g_vertex_buffer_data[3*i+1] = massPositionFloat(cubeVertex(triangleVertex(i))+k,1);       ``````````````````
//        g_vertex_buffer_data[3*i+2] = massPositionFloat(cubeVertex(triangleVertex(i))+k,2);
//        //std::cout<<massPosition(cubeVertex(triangleVertex(i)),0)<<massPosition(cubeVertex(triangleVertex(i)),1)<<massPosition(cubeVertex(triangleVertex(i)),2)<<std::endl;
//    }
//


    GLfloat g_color_buffer_data_template[] = {
            0.583f,  0.771f,  0.014f,
            0.609f,  0.115f,  0.436f,
            0.327f,  0.483f,  0.844f,
            0.822f,  0.569f,  0.201f,
            0.435f,  0.602f,  0.223f,
            0.310f,  0.747f,  0.185f,
            0.597f,  0.770f,  0.761f,
            0.559f,  0.436f,  0.730f,
            0.359f,  0.583f,  0.152f,
            0.483f,  0.596f,  0.789f,
            0.559f,  0.861f,  0.639f,
            0.195f,  0.548f,  0.859f,
            0.014f,  0.184f,  0.576f,
            0.771f,  0.328f,  0.970f,
            0.406f,  0.615f,  0.116f,
            0.676f,  0.977f,  0.133f,
            0.971f,  0.572f,  0.833f,
            0.140f,  0.616f,  0.489f,
            0.997f,  0.513f,  0.064f,
            0.945f,  0.719f,  0.592f,
            0.543f,  0.021f,  0.978f,
            0.279f,  0.317f,  0.505f,
            0.167f,  0.620f,  0.077f,
            0.347f,  0.857f,  0.137f,
            0.055f,  0.953f,  0.042f,
            0.714f,  0.505f,  0.345f,
            0.783f,  0.290f,  0.734f,
            0.722f,  0.645f,  0.174f,
            0.302f,  0.455f,  0.848f,
            0.225f,  0.587f,  0.040f,
            0.517f,  0.713f,  0.338f,
            0.053f,  0.959f,  0.120f,
            0.393f,  0.621f,  0.362f,
            0.673f,  0.211f,  0.457f,
            0.820f,  0.883f,  0.371f,
            0.982f,  0.099f,  0.879f
    };

    Eigen::ArrayXf g_color_buffer_data(num_cubes*num_vertices_per_cube*3,1);
    for (int m = 0; m < g_color_buffer_data.size(); ++m) {
        g_color_buffer_data.data()[m] = g_color_buffer_data_template[m%(num_vertices_per_cube*3)];
    }





    GLuint vertexbuffer;
    glGenBuffers(1, &vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
//    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);


    GLuint colorbuffer;
    glGenBuffers(1, &colorbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, g_color_buffer_data.size()* sizeof(float), &g_color_buffer_data.data()[0], GL_STATIC_DRAW);


///////////////////////////////////////buffers for the floor////////////////////////////////////////////////////////////
    GLuint VertexArrayID2;
    glGenVertexArrays(1, &VertexArrayID2);
    glBindVertexArray(VertexArrayID2);


    int floor_length = 1000;

    int num_squares_x = 55;
    int num_squares_y = 55;
    int num_squares = num_squares_x*num_squares_y;

    Eigen::Array<float,Eigen::Dynamic,3,Eigen::RowMajor> vertexbuffer_floor_data_square(6,3);
    vertexbuffer_floor_data_square<<
            0,0,0,
            1,0,0,
            0,1,0,
            1,1,0,
            1,0,0,
            0,1,0;

    Eigen::Array<float,Eigen::Dynamic,3,Eigen::RowMajor> vertexbuffer_floor_data(num_squares*6,3);





    Eigen::Array<float,Eigen::Dynamic,3,Eigen::RowMajor> colorbuffer_floor_data_square(6,3);
    colorbuffer_floor_data_square<<
            255.0 / 256, 77.0 / 256, 77.0 / 256,
            255.0 / 256, 10.0 / 256, 77.0 / 256,
            255.0 / 256, 77.0 / 256, 77.0 / 256,
            255.0 / 256, 40.0 / 256, 77.0 / 256,
            55.0 / 256, 77.0 / 256, 77.0 / 256,
            255.0 / 256, 77.0 / 256, 77.0 / 256;

    Eigen::Array<float,Eigen::Dynamic,3,Eigen::RowMajor> colorbuffer_floor_data(num_squares*6,3);



#pragma omp parallel for collapse(2)
    for (int i= 0; i <num_squares_x; ++i) {
        for (int j = 0; j <num_squares_y; ++j) {
            int square_k = num_squares_y*i+j;
            vertexbuffer_floor_data.block(square_k*6,0,6,3) = vertexbuffer_floor_data_square;
            vertexbuffer_floor_data.block(square_k*6,0,6,1) += i- double(num_squares_x)/2;
            vertexbuffer_floor_data.block(square_k*6,1,6,1)+= j- double(num_squares_y)/2;

            colorbuffer_floor_data.block(square_k*6,0,6,3) = colorbuffer_floor_data_square;
        }
    }


    GLuint vertexbuffer_floor;
    glGenBuffers(1, &vertexbuffer_floor);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer_floor);
    glBufferData(GL_ARRAY_BUFFER, vertexbuffer_floor_data.size()*sizeof(float), &vertexbuffer_floor_data.data()[0], GL_STATIC_DRAW);

    GLuint colorbuffer_floor;
    glGenBuffers(1, &colorbuffer_floor);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer_floor);
    glBufferData(GL_ARRAY_BUFFER, colorbuffer_floor_data.size()*sizeof(float), &colorbuffer_floor_data.data()[0], GL_STATIC_DRAW);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//    // Generate a buffer for the indices as well
//    GLuint elementbuffer;
//    glGenBuffers(1, &elementbuffer);
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
//    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned short), &indices[0] , GL_STATIC_DRAW);

// fps counter : http://www.opengl-tutorial.org/miscellaneous/an-fps-counter/


    int num_frames = position_history.rows()/num_masses;

    double lastTime = glfwGetTime();
    int nbFrames = 0;
    int k_frames = 0;
    double dt_render = 1.0/60.0;

    do{

        processInput(window);
        view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
        MVP  = Projection * view * Model;


        // Measure speed
        double currentTime = glfwGetTime();
        if(currentTime - lastTime >= dt_render){
            if ( currentTime - lastTime >= 1.0 ){ // If last prinf() was more than 1 sec ago
                // printf and reset timer
                printf("%f ms/frame\n", 1000.0/double(nbFrames));
                nbFrames = 0;
                lastTime += 1.0;
            }



            for (int i = 0; i < num_masses; ++i) {
                for (int j = 0; j < 3; ++j) {
                    massPositionFloat(i,j) = (float)position_history(k_frames*num_masses+i,j);
                }
            }


            for (int k = 0; k < num_cubes; ++k) {
                int cube_position = springtoMass(28*k,0);
                //std::cout<<cube_position<<std::endl;
                //draw 36 vertices
                for (int i=0;i<num_vertices_per_cube;i++){

                    g_vertex_buffer_data[num_vertices_per_cube*k*3+3*i+0] = massPositionFloat(cubeVertex(triangleVertex(i))+cube_position,0);
                    g_vertex_buffer_data[num_vertices_per_cube*k*3+3*i+1] = massPositionFloat(cubeVertex(triangleVertex(i))+cube_position,1);
                    g_vertex_buffer_data[num_vertices_per_cube*k*3+3*i+2] = massPositionFloat(cubeVertex(triangleVertex(i))+cube_position,2);
                    //std::cout<<massPosition(cubeVertex(triangleVertex(i)),0)<<massPosition(cubeVertex(triangleVertex(i)),1)<<massPosition(cubeVertex(triangleVertex(i)),2)<<std::endl;
                }
            }


            glBindVertexArray(VertexArrayID);
            glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
            glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);




            // Clear the screen
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // Use our shader
            glUseProgram(programID);

            // Send our transformation to the currently bound shader,
            // in the "MVP" uniform
            glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

            // 1rst attribute buffer : vertices
            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
            glVertexAttribPointer(
                    0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
                    3,                  // size
                    GL_FLOAT,           // type
                    GL_FALSE,           // normalized?
                    0,                  // stride
                    (void*)0            // array buffer offset
            );

            // 2nd attribute buffer : colors
            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
            glVertexAttribPointer(
                    1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
                    3,                                // size
                    GL_FLOAT,                         // type
                    GL_FALSE,                         // normalized?
                    0,                                // stride
                    (void*)0                          // array buffer offset
            );

            // Draw the triangle !
            glDrawArrays(GL_TRIANGLES, 0, 12*3*num_cubes); // 12*3 indices starting at 0 -> 12 triangles

            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);

///////////////////////////attribute buffer for the floor/////////////////////////////////////////////////////////
            glBindVertexArray(VertexArrayID);
            glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer_floor);
            glVertexAttribPointer(
                    0,                                // attribute. No particular reason for 1, but must match the layout in the shader.
                    3,                                // size
                    GL_FLOAT,                         // type
                    GL_FALSE,                         // normalized?
                    0,                                // stride
                    (void*)0                          // array buffer offset
            );

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorbuffer_floor);
            glVertexAttribPointer(
                    1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
                    3,                                // size
                    GL_FLOAT,                         // type
                    GL_FALSE,                         // normalized?
                    0,                                // stride
                    (void*)0                          // array buffer offset
            );


            glDrawArrays(GL_TRIANGLES, 0, 3*2*num_squares); // 12*3 indices starting at 0 -> 12 triangles

            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            // Swap buffers
            glfwSwapBuffers(window);
            glfwPollEvents();

            nbFrames++;
            k_frames++;
            if(k_frames>=num_frames){
                k_frames = 0;
            }
        }




    } // Check if the ESC key was pressed or the window was closed
    while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           glfwWindowShouldClose(window) == 0 );

    // Cleanup VBO and shader
    glDeleteBuffers(1, &vertexbuffer);
    glDeleteBuffers(1, &colorbuffer);
    glDeleteBuffers(1, &vertexbuffer_floor);
    glDeleteBuffers(1, &colorbuffer_floor);
    glDeleteProgram(programID);
    glDeleteVertexArrays(1, &VertexArrayID);

    // Close OpenGL window and terminate GLFW
    glfwTerminate();

    return 0;
}



double simulate(
        int num_frames,
        int skip_frames,
        ArrayXdRowMajor& springCoefa,
        ArrayXdRowMajor& springCoefb,
        Eigen::ArrayXd& l0,
        ArrayX3dRowMajor& massPosition,
        ArrayX2dRowMajor& springtoMass,
        ArrayX3dRowMajor& massVelocity,
        ArrayX3dRowMajor& massAcceleration,
        ArrayX3dRowMajor& position_history
        ){

    double time_stamp = 0.000;
    ArrayX3dRowMajor massForces = ArrayX3dRowMajor::Zero(num_masses,3);
    for (int i = 0; i <num_frames; ++i) {
        for (int j = 0; j <skip_frames; ++j) {
            massForces.setZero();

            Eigen::ArrayXd l0t = Setl0t(
                    time_stamp,
                    springCoefa,
                    springCoefb,
                    l0);

            applyGravity(massForces);
            applySpringForces(l0t,massPosition,massForces,springtoMass);
            applyFrictionandGroundForces(massPosition,massForces,massVelocity);
            massAcceleration = massForces/Mass;
            massVelocity += massAcceleration*TimeStep;
            massVelocity *= DampingCoef;
            massPosition += massVelocity*TimeStep;

            time_stamp += TimeStep;
        }
        position_history.block(i*num_masses,0,num_masses,3) = massPosition;
    }
    double fitness = 0;
    for (int i=0;i<c+1;i++){
        for (int j=(a+1)*b;j<(a+1)*(b+1);j++)
        fitness += massPosition(i*(a+1)*(b+1)+j,1);
    }
    fitness /= (a+1)*(c+1);
    return fitness;
}

template<typename T>
Eigen::ArrayXi orderFitness(const Eigen::DenseBase<T>& values) {
    Eigen::ArrayXi indices(values.size());
    for (int i = 0; i <values.size(); ++i) {
        indices(i) = i;
    }
    std::sort(
            indices.begin(), indices.end(),
            [&](size_t a, size_t b) { return values(a) > values(b); }
    );
    return indices;
}



void savefitness(
        int generation_i,
        double best_fitness,
        std::string file_path){
    std::ofstream myfile;
//    if (generation_i==0)
//        myfile.open(file_path);
//    else
    myfile.open(file_path,std::ios_base::app);
    myfile<<generation_i<<"\t"<<best_fitness<<std::endl;
//    myfile<<best_fitness<<std::endl;
    myfile.close();
}



void savefitnesshistory(
        int generation_i,
        Eigen::ArrayXd best_fitness_arr){
    std::ofstream myfile;
    if (generation_i==0)
        myfile.open("fitness_history.txt");
    else
        myfile.open("fitness_history.txt",std::ios_base::app);
    myfile<<best_fitness_arr<<std::endl;
    myfile.close();
}
int main() {



    Eigen::ArrayXi cubeVertex = initCubeVertex();
    Eigen::ArrayXi triangleVertex = initTriangleVertex();
    ArrayX3dRowMajor baseMassPosition = initBaseMassPosition();
    ArrayX3dRowMajor massVelocity = initMassVelocity();
    ArrayX3dRowMajor massAcceleration = ArrayX3dRowMajor::Zero(num_masses, 3);
    ArrayX2dRowMajor baseSpringtoMass = initBaseSpringtoMass();
    ArrayX2dRowMajor springtoMass = initSpringtoMass(
            baseSpringtoMass,
            baseMassPosition);
    Eigen::ArrayXd l0 = initL0();
    ArrayX3dRowMajor massPosition = initMassPosition();


////////////////////// evolve ////////////////////////////////////////////////////
//    int num_frames = 500;
//    int skip_frames = 32;
//    int num_evaluations = 2048*64;
//    int population_size = 64*4;
//    int selection_pressure = 64*2;
//
////    int num_evaluations = 64*10;
////    int population_size = 16;
////    int selection_pressure = 8;
//
//    int num_generation = (num_evaluations- population_size)/selection_pressure+1;
//
//    Eigen::ArrayXd best_fitness_arr(population_size);
//
//    ArrayXXdRowMajor springCoefa_arr = initSpringCoefa(population_size,num_cubes);
//    ArrayXXdRowMajor springCoefb_arr = initSpringCoefb(population_size,num_cubes);
//
//
//    std::random_device rd;  //Will be used to obtain a seed for the random number engine
//    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//    std::uniform_real_distribution<> dis(0.0, 1.0);
//    std::uniform_int_distribution<> dis_selection(0,selection_pressure-1);
//    std::uniform_int_distribution<> dis_cubes(0,num_cubes-1);
//    std::uniform_real_distribution<> dis_coefa(0.0, CoefaRange);
//    std::uniform_real_distribution<> dis_coefb(0.0, CoefbRange);
//
//    std::ostringstream ss;
//    ss<< dis(rd);
//    std::string random_string(ss.str());
//    std::string file_name_diversity = "./data/diversity_"+random_string+".txt";
//    std::string file_name_fitness = "./data/fitness_"+random_string+".txt";
//
//
//    int best_solution = 0;
//    double best_fitness = -999;
//    ArrayX3dRowMajor best_position_history= ArrayX3dRowMajor::Zero(num_frames*num_masses,3);
//
//
//    //the first generation
//#pragma omp parallel for
//    for (int i = 0; i < population_size; i++) {
//        //std::cout<<i<<std::endl;
//        //
//        springCoefa_arr.block(i, 0, 1, num_cubes)<<-0.0787365553522187,0.0397860520993408,0.0285811549930746,0.0148356401616873,0.0225603931014757,0.00942131270161891,0.0190325777693803,0.0233728654418946,0.0583954952556619,-0.0483040296883807,0.0242124729536801,0.00489416932481162,0.0665881243375549,0.0251908738910303,0.0741066694232201,0.0670986463979148,-0.0345643036600967,-0.0259200940029323,-0.0768512000128865,-0.066696290442113,0.0769517727367604,0.0440311868354987,0.00310712193810026,0.078899163230741,0.0788787391934348,0.0771249259249889,0.0117853437113109,0.0767958074048142,0.0437036912719271,0.0292698902400536,0.0590928833160164,0.0449671359942142,-0.0638690479769693,0.0677258177135295,0.0741325108270997,0.0139300205356702,-0.0556104941552554,0.078835353850776,0.00860470357830995,-0.0199256604630154,0.0330242166540873,0.0383049539282475,0.0671731600105638,0.0100578547688471,-0.0763951291127108,0.0609859394193561,0.0691503421495137,0.0247812721674833,0.0637547364010265,-0.0723889178002202,-0.0519388390946849,0.0624782460846371,-0.0431796303219998,0.00900752563448974,-0.0400732179545207,0.0376800872188434,0.021020399800878,0.0117542496590411,0.0336049248422409,0.0625928757258658,0.0278127577633072,0.0702474944993143,0.0435802303085012,-0.0415860856145556,0.0662081225558886,-0.056616993684609,0.0702688369854674,0.0669085326922944,-0.0577816398338329,0.00498168028038696,0.064275964950129,0.0082106890540437,-0.0675601910741814,-0.0297229902212149,0.0547141934895084,-0.0238513961172902,0.0504024790089589,-0.0269918353608771,0.0652930952540101,0.0246078541616946,0.0783054675153948,0.0217986215452878,-0.0539794790623614,-0.0749792838073239,0.0455257733005685,0.0198429829005026,0.0347960266260412,-0.0795648236570716,-0.060654577007822,-0.0186328556102854,0.079048270806227,0.0673910252187208,-0.0283853611854768,0.0768474262419641,4.22997400361602e-05,-0.0458803351064587,0.0660115074301192,-0.00968886334900225,0.047408684808348,-0.0717701324782195,0.035807761529371,-0.0174382546997807,0.0189957893067148,0.0482475704551896,0.0328387551535101,0.0190584340176443,-0.0556038256621006,0.00324123408796323,0.0127268029229173,0.0472632276808765,-0.052150911824848,0.0281599258042537,0.0174305656260953,-0.0261303908872094,0.035257932212356,-0.017043661147842,0.0726981328223454,-0.0570304115009636,0.0253911636689594,0.0197565364348446,0.016967082249871,0.0102833213978636,-0.0188242315961161,0.0559513717032743,-0.0549317130329701,0.0612180681439201,0.0193526401982814,-0.0689202056773567,0.0065515928256988,0.00290657720210824,-0.0606903381555762,-0.0726630338247228,-0.0690636828490829,0.0193897279593573,0.0555845366304668,0.0437750723044272,0.0468450873002075,0.0394696528074704,-0.0329836936821154,0.0769758056090101,-0.0703300195887359,-0.00513460550696338,0.0571396772193301,0.0142195497380376,0.0164515214097471,0.0308472113656125,0.0694353171229428,0.0780693864068339,0.0725435366842717,0.00753803458415812,-0.0130345523790617,0.0691651284023265,-0.0578139640287561,-0.0610438296855631,-0.0769839821741842,-0.00349377188994259,-0.0491218169634798,0.0680368567239194,0.0104728014513025,0.0277090399468825,0.0643284656034449,0.0148744277499296,0.0319480187967178,0.0469556046268861,0.0225398461253102,0.0540060857562377,-0.0767194037869197,-0.0488786590699473,0.0358852142383299,0.0235482812711195,0.00388935815724048,-0.0344991294827774,0.0072645609730069,0.0182211471806379,-0.0330087505062151,0.066253237997802,-0.007020162309995,0.0507489917067276,0.0579286789791326,-0.0554916582887488,-0.0443672990633023,-0.0351058734744349,0.0638490229759595,0.062662030843652,-0.016149703159998,0.0502813330899371,0.054324965017999,0.0552459516260987,0.0749330744627788,0.075759350480195,0.058129147798243,0.0282569440958355,-0.0268196486434059,-0.022496843981788,0.0128301607010894,0.00652780320799341,0.00679419323244019,0.0513311464426714,-0.0737902242661408,0.00482908282653852,-0.028651634505322,-0.00741968561309375,-0.0113375247695192,-0.0309082336495203,-0.0691985385069617,0.0356537247242656,-0.0510489410772216,0.00378129925754912,-0.0162516089231947,-0.0731202621725948,0.0282896409688004,0.019381092013503,0.00626741023776652,0.0410187031411662,0.0771998289214446,0.0356241612674781,0.00990672205602387,-0.0794443978366649,0.0199252965021531,-0.0576569127932456,0.00239666443429732,0.000663004027198092,-0.0266708465417246,0.0555770157908914,0.0329926150670633,0.032585694171761,-0.017895181075621,0.070577879464054,0.0374795822822829,-0.0116854053417618,-0.00815244164697474,0.00814443644515446,0.0608949091196502,0.0605100336580118,0.05723620287014,0.0716963706126885,0.0161637583077716,-0.0738127382815875,-0.00452233020426814,0.0799121494590827,-0.0669330003796765,0.0132565554729351,0.0192932413980799,-0.0151591360267061,0.0587106682500125,0.0164930703195245,-0.0595349748337339,0.0651585560967953,-0.0119821358900434,0.0313351379853371,0.0102307653847294,0.0735335526123895,-0.0783080907903184,0.00670429681703182,0.0459915443351453,0.029551836526744,0.01614561294026,-0.0519036367404757,0.0165703120532307,0.0635226350869781,0.0793850834321845,-0.0715821296682498,0.00138656943146689,0.0789834328713299,0.068927903989762,-0.0416776768684749,0.0791554930915526,0.00509166229753367,-0.0354904150755566,0.0487799919627843,0.0385185939230934,0.0429319663492875,-0.0717527819386464,-0.0557029469943153,0.0424174485180608,0.0275694849097959,0.0283013731032796,0.0628824737588328,0.065324034305884,0.00756486963832979,-0.0580133232371944,0.000926739197192052,0.040198942275384,0.0671719744676738,-0.0315740670038266,0.00459441494531404,0.0406146016064966,-0.0549752013268765,0.0632378235381273,0.0168862972014782,0.0632471188079785,0.0675372461255414,0.00629654168699845,-0.0716245620658736,0.0733686844228621,-0.0731565964935145,0.0207637777727245,-0.063657464228411,0.0332426104476874,0.0721893606780537,-0.059819667255422,0.0611593644201538,0.0371006349274426,0.008756611388529,-0.0662300293828501,0.0403580811924986,0.0173801001790549,-0.0183408645998411,0.0773245660948216,0.0381073966427275,-0.0354915834756063,0.0505520611814888,0.0779140615453543,0.0467182230794422,0.00691535268300928,-0.00752095301054463,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.072972533086243,-0.0563313566038065,0.00493563063820985,0.0168609884832338,-0.0755998318621888,-0.0764228115959199,0.0651290465822113,0.0677687595727708,-0.0255279210701249,0.0154125595537073,0.0310159908752032,0.0324546150026772,0.0672664034017367,-0.014407729028914,0.0427632872586899,-0.0203874216696189,-0.0329939846475581,0.0582008473473605,0.049530052081463,0.0420176874102258,-0.0526738212130376,0.0566280047414426,-0.0299665013467737,0.0202939150204388,-0.00878208395502627,0.0613053419223006,0.0651859232358668,0.0415810448984855,0.0684167652523223,0.0225791594861915,-0.0764551426453773,-0.00469409365425543,0.0657973458551789,-0.0340542982304722,-0.0545746671290019,-0.0363075187598856,-0.0704504885480043,0.0251474474115052,0.0239696648456015,0.00888312433328625,-0.0184900807302865,-0.026193907384851,0.0346668059912821,-0.0502069033543611,0.0727939190186036,-0.0652426554193919,-0.073405683857112,0.058832445019454,0.0433435391706107,-0.0276349627727805,-0.0627231468924895,0.0216381766375332,0.0484023522030692,0.0511634874649726,0.0491580558268187,-0.0303264040268615,0.0372968669781912,0.0201023868931081,0.0588998787006828,0.0257136321560077,-0.0356257363947228,0.0141284905145072,-0.0589804615727535,-0.0498283906140497,-0.0516095621751666,-0.0335551287017554,-0.00613590937393527,-0.0420600507231709,0.0715923187097499,-0.0621662446028396,-0.0562829418369955,-0.0100930846902044,-0.00836015198769055,0.00148987960139749,0.00289533462510226,0.0398096824344578,0.0162472241820056,0.00948965076799023,-0.0518638096618763,0.0600174140898688,0.0472304060902588,-0.0355273761393164,0.0202883375903072,0.0234392278378081,0.0422216261339447,-0.0518900550770993,0.0370514884947964,0.0149342752690121,0.0455472842443489,0.0482517623753528,0.0192839286938701,0.05911724067252,-0.00100677106576356,-0.0409349269237998,0.0659956001425141,0.0674249651969527,-0.0747426529204206,0.0431735729416369,0.0108531262403601,0.00681786573268896,4.67374548533606e-05,0.0628869602749529,-0.0176661403186927,0.0786480859660767,0.0524744645035985,-0.0418882346720846,0.00538064945739723,0.0144900607737209,0.0412227917798612,0.0745418645046977,0.0701532715513153,0.000597512655533077,-0.0725604402984308,0.0226924941421917,0.0577535415337205,0.0666303270992964,0.0287497000157599,-0.0541356280325612,0.0614123152109852,-0.0198186871129175,0.0215594750847957,0.0270599850426697,0.0423853336672181,-0.0427899707680521,0.0499074528313328,0.0793047787834378,0.0473748562227119,0.0524015135663986,-0.0379093903107147,0.0793937288558049;
//        springCoefb_arr.block(i, 0, 1, num_cubes)<<-2.43192959884493,6.08514030121006,4.04035678062355,4.0939778826623,0.885759964104128,-4.40632201143234,-3.84532036998163,-5.33669165449445,-0.201902543387136,-5.74421405700816,1.84616840882136,3.98087731087271,1.11369296464456,0.934977743915378,2.02614908115595,0.959073349058476,5.15053093739066,3.006446367785,-1.6732576622727,-1.23768408693475,1.27684178065925,6.04271727318241,3.77882575946868,-6.08085982355576,0.652931759876538,-5.48543869292053,1.81651248891278,5.96869830004477,1.31094790673134,3.23199520554859,0.163817214486977,5.16220350788641,-2.7309598098749,5.48964219698341,5.39277565112383,3.98484668574114,1.44712120543646,5.41086102056708,5.39994485533409,0.946468807561242,4.0832679471894,0.279088417230802,-4.91215907281787,4.6082018532571,0.92838679515286,-1.64664393517429,3.42329263641773,5.70322596816668,-2.65748894621183,1.29342632652648,0.568687125340165,-2.43705856028738,-1.92000280751018,-1.8801132043928,-0.517201690556528,1.68458566290145,4.00263656362903,2.50000791231667,0.452871199207279,-4.12575887624951,1.3802421085177,2.0898833298572,1.09336188408456,4.51445899470899,0.99940722926591,-2.21682713707008,1.5502626395456,5.9937160473333,-3.08915132235467,5.79470203052607,5.86985869733578,1.12879047483547,2.60477802497352,0.454742152576666,0.517470415219057,-5.5002405095043,-0.520598258586578,-5.43747913056799,-2.17154290659146,-4.77267751366817,0.0886738094516904,2.80624253330262,6.11985611039204,-5.33693277264827,-6.2580654340777,1.58892567487088,1.83628608125866,-5.35466976461089,2.3785831099318,-3.87135237064131,3.14285901158726,4.07157566778108,4.50171615921588,4.25724858523384,-4.71109470247867,-0.420077554397874,2.01939375860174,3.12235313706693,3.04526415542555,5.21342763624706,-3.62827223889698,-0.737268266315717,1.96095282746985,5.25969098607654,6.00065908626095,3.66813714062273,6.04263567657224,-0.803124378177302,1.71191742881433,0.288838383269688,0.707383314006201,1.04678416576555,-2.49262183233593,0.54405421854657,4.53851627650281,-2.46750206641364,2.16747008949913,-0.845862873766654,2.96766701262683,1.696825034699,1.26885943120368,1.67349275491913,-1.65457166366446,-0.215499085192086,-0.443091929617129,-0.0824811602914577,3.72336076311685,-4.70688337686707,4.86531505313461,1.28618349288566,-5.77664094062,-0.588400267973177,0.716840377031068,4.52819929759633,-1.61189448189663,0.434314257440346,4.15175172698769,2.17904391568207,5.91437508511471,0.300796283741831,2.01864875880472,0.338573193269244,5.38597896424107,5.14448955286199,3.82100555758673,4.50086791031566,3.6169715871486,3.34167081581541,5.25594889042253,3.6242449893384,2.869005561767,5.92319047857269,4.68127868871957,-3.14173898982101,-2.97035686401395,-0.230860614223697,-0.657414860877198,0.394065895226855,1.19275490160846,5.56181821803563,-2.91342920436725,1.82579686459939,-1.63544369410512,1.79662060059091,3.02465413262939,-1.00109808835509,4.72917326240044,4.64586125535466,2.46719717703016,5.95293752900116,-0.728069876461334,5.10586141337541,2.35273521470591,1.07370842850819,-1.93135481169837,4.37060868490869,-5.19487679371209,0.391702183144271,-1.21653276126522,-6.13711309553413,-3.61258757467371,-4.63071240534989,1.85101678249404,0.308870603002905,-1.4892661951709,5.49822140170814,0.837830493970494,5.16937915605619,5.61545234316673,3.98285231687726,4.89567875026041,-6.23992361757923,-5.34797394185981,-2.9636930323549,1.39715120233631,-2.12644323785466,3.40947531887144,1.80457371889052,4.09084134539069,-4.52002808382579,1.6512736377946,5.63214168372345,3.85473970142402,4.58615899606649,0.42266490637997,-4.35980031027435,0.0279554131341941,1.51097331266788,0.96973753156521,5.09460785186898,1.65704541713375,3.6403351568915,1.05442214408828,0.661927881806177,4.70902626508569,-1.02537075450348,0.192632458299389,-2.51630261939429,1.84969643455691,3.83759113021315,1.5746427306926,5.27170062187801,-4.70822915029589,2.50985398883279,2.6840672558054,-6.16522383252989,-5.89977445487354,3.05056480308843,2.21451166528256,4.47425209051715,5.5492140964325,3.75304063429704,3.82320856838893,3.12076859200485,2.05601442451186,-2.03731173108278,5.0441534817305,-4.19921536235395,5.75684678158511,-0.269294192555964,-5.38779271633664,5.43889699630295,-2.91214424151613,-4.92389726981756,3.36810651670152,-4.48630318228211,0.33391717567896,-4.53575025726799,1.0605535709224,4.47792557723166,-3.00405579583632,3.00532537434112,-0.591793793133674,2.71423009991725,5.14505028459612,-2.34202100170558,-1.31113837281043,5.52846103557426,-4.83914922729897,3.79489135958476,3.12888332563001,-5.57312033671814,1.83308250263569,3.19973161421444,3.83083346113838,-5.96668672335076,3.33493105350654,2.5918017370172,-3.88271688570471,3.38475512344542,4.01203987647953,2.89624094641009,-0.173391873853497,-0.589636697054897,-1.62803646600724,-1.48343160398213,0.848696330955331,4.98906590967172,4.47871139269655,1.19791773758934,-4.04150873496511,3.90776009537591,3.43991137256052,2.61626586932033,1.49626440893562,2.84893529190212,4.7158110969248,-3.53724639953568,4.50491096743946,1.49275078504112,-5.77939771606788,4.17997546079867,4.1568390815018,-6.27213179494185,-5.02828562505671,-3.57958680001224,3.34456559070511,3.46028085870612,0.799359763371303,1.47706369856494,2.56591515626386,3.63946459053563,5.86683490040684,5.73396412750752,4.27641316775727,4.53509332871182,5.61569700322845,2.9393886033297,-2.31855975768557,-5.79668449207325,-5.13523515189644,4.72865749192539,-1.49676813258366,4.1860763159865,3.37062530542857,1.21079109463131,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,5.12690620844648,3.61694724017608,3.05138930367418,-4.01275011252449,-2.72372981137355,-0.0180124155448347,4.12303341309968,-3.74873780023526,-0.887138111904024,-5.41990213502721,-5.9774984601979,0.552717896425435,0.921612258771997,3.17838284590814,1.86797359438216,-1.47430087116196,1.35138812294626,-3.90908982234295,0.382622065522161,3.58663702089682,-1.89477216410266,0.842265657489956,5.32264774649661,4.07279963344888,5.93452553491607,4.58997732369303,0.786062578409914,0.325200986907274,-5.9917002755706,0.588543628002016,-3.55057119553441,1.28408373976525,-1.25341481552174,-2.18636580089857,1.02549154472145,-2.97758908615075,-3.02601767758674,-0.79860588073586,-1.31513570862061,-5.10264541624672,-4.04774808632133,4.46998918241493,1.57982980509186,5.8374249517381,0.290319230405863,1.86919610233567,0.185816115052586,2.65571266392906,3.43118334743,1.54324804109866,-0.19298985557722,-5.47214491078457,2.00854924420448,2.55367598272938,2.67066286720342,-4.40899082454341,-4.04429312445419,5.6046481394415,-1.30713353515645,-3.75280819417312,5.43198984070525,0.316535709131635,3.81446074559214,-2.10461017481648,5.52229986841057,-1.44323290968642,1.20098593903276,-3.78690300917618,4.04134640957772,6.16903543626382,-0.65473740484626,-5.61395256815843,4.35583941282708,5.25665177966897,5.8310232749945,1.77763132207111,0.842662682004633,-0.266345809952913,1.76766768752563,2.41652541542875,4.29474875067796,5.02212468610824,-3.23113956818491,-4.46801644912111,3.88677049527466,-5.81971337358793,2.74959988627578,-3.37151591823711,2.76592180759093,-3.16280122963266,2.80490867323494,2.55004818904424,0.689441449147486,-0.0351284046447834,6.06388826542954,3.55993961157682,0.535486416593954,3.05170925752847,-1.00278238102184,5.46116454437677,5.81344950372511,-1.74189744386769,-0.332778030647898,-4.03931177144349,4.25502161714161,4.25874652459808,1.4159745097112,3.55494616590075,4.0525786915951,5.09901468640174,4.31990406550899,1.04140371424045,5.82363717143333,-2.31027463831456,-0.906303222026946,-0.757257994087879,-6.03794887918608,0.696102823104212,1.00660345861309,4.57369870271486,5.17365517040888,0.38104645055342,3.68186372776027,3.74661889802106,3.58529704079631,4.04925640652263,4.59815072532965,1.95431786572231,-5.33829054200515,1.63610028456171;
//
//
//        ArrayXdRowMajor springCoefa_i = springCoefa_arr.block(i, 0, 1, num_cubes);
//        ArrayXdRowMajor springCoefb_i = springCoefb_arr.block(i, 0, 1, num_cubes);
//
//        ArrayX3dRowMajor massPosition_i = ArrayX3dRowMajor(massPosition);
//        ArrayX3dRowMajor massVelocity_i = ArrayX3dRowMajor::Zero(num_masses, 3);
//        ArrayX3dRowMajor massAcceleration_i = ArrayX3dRowMajor::Zero(num_masses, 3);
//
//        ArrayX3dRowMajor position_history_i = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);
//
//        best_fitness_arr(i) = simulate(
//                num_frames,
//                skip_frames,
//                springCoefa_i,
//                springCoefb_i,
//                l0,
//                massPosition_i,
//                springtoMass,
//                massVelocity_i,
//                massAcceleration_i,
//                position_history_i);
//        if(i==0){
//            best_position_history =position_history_i;
//            best_solution = i;
//        }
//    }
//
//
//
//    int best_i = orderFitness(best_fitness_arr)(0);
//
//    Eigen::IOFormat print_formating(Eigen::FullPrecision,Eigen::DontAlignCols,",",",","","","",";");
//
//    std::cout<<"springCoefa<<"<<springCoefa_arr.row(best_i).transpose().format(print_formating)<<std::endl;
//    std::cout<<"springCoefb<<"<<springCoefb_arr.row(best_i).transpose().format(print_formating)<<std::endl<<std::endl;
//
//
//    best_fitness = best_fitness_arr(best_i);
//
//    int  num_evaluated= population_size;
//    savefitness(num_evaluated,best_fitness,file_name_fitness);
//
//
//
////    #pragma omp parallel for
//    for (int i_generation = 1; i_generation < num_generation; i_generation++) {
//        Eigen::ArrayXi fitness_indices(population_size);
//
//        fitness_indices = orderFitness(best_fitness_arr);
//        std::cout << "generation" <<i_generation<< std::endl;
//
//        #pragma omp parallel for
//        for (int j = 0; j < selection_pressure; ++j) {
//            // crossover
//            int parent1 = fitness_indices(j);
//            int parent2;
//            do {
//                parent2 = fitness_indices(dis_selection(gen));
//            } while (parent1 == parent2);
//            int st, ed;
//            do {
//                st = dis_cubes(gen);
//                ed = dis_cubes(gen);
//            } while (st > ed);
//
//            ArrayXdRowMajor springCoefa_i = springCoefa_arr.block(parent1, 0, 1, num_cubes);
//            ArrayXdRowMajor springCoefb_i = springCoefb_arr.block(parent1, 0, 1, num_cubes);
//
//            springCoefa_i.block(0, st, 1, ed - st + 1) = springCoefa_arr.block(parent2, st, 1, ed - st + 1);
//            springCoefb_i.block(0, st, 1, ed - st + 1) = springCoefb_arr.block(parent2, st, 1, ed - st + 1);
//
//
//            // mutation
//            int mutation_point = dis_cubes(gen);
//
//            springCoefa_i(mutation_point) = dis_coefa(gen);
//            springCoefb_i(mutation_point) = dis_coefb(gen);
//
//
//            ArrayX3dRowMajor massPosition_i = ArrayX3dRowMajor(massPosition);
//            ArrayX3dRowMajor massVelocity_i = ArrayX3dRowMajor::Zero(num_masses, 3);
//            ArrayX3dRowMajor massAcceleration_i = ArrayX3dRowMajor::Zero(num_masses, 3);
//
//
//            ArrayX3dRowMajor position_history_i = ArrayX3dRowMajor::Zero(num_frames * num_masses, 3);
//
//            double fitness_child = simulate(
//                    num_frames,
//                    skip_frames,
//                    springCoefa_i,
//                    springCoefb_i,
//                    l0,
//                    massPosition_i,
//                    springtoMass,
//                    massVelocity_i,
//                    massAcceleration_i,
//                    position_history_i);
//
//
//            if (fitness_child > best_fitness_arr(parent1)) {
//                //replace this parent
//                springCoefa_arr.block(parent1, 0, 1, num_cubes) = springCoefa_i;
//                springCoefb_arr.block(parent1, 0, 1, num_cubes) = springCoefb_i;
//                best_fitness_arr(parent1) = fitness_child;
//
//                if (fitness_child > best_fitness) {
//                    best_i = parent1;
//                    best_position_history = position_history_i;
//                    best_fitness = fitness_child;
//                    std::cout<<"best = "<<best_fitness<<std::endl;
//                    std::cout<<"springCoefa<<"<<springCoefa_arr.row(best_i).transpose().format(print_formating)<<std::endl;
//                    std::cout<<"springCoefb<<"<<springCoefb_arr.row(best_i).transpose().format(print_formating)<<std::endl<<std::endl;
//                }
//            } else {
//                int parent_for_replacement = fitness_indices(selection_pressure + j);
//                springCoefa_arr.block(parent_for_replacement, 0, 1, num_cubes) = springCoefa_i;
//                springCoefb_arr.block(parent_for_replacement, 0, 1, num_cubes) = springCoefb_i;
//                best_fitness_arr(parent_for_replacement) = fitness_child;
//                // put it somewhere else
//            }
//        }
//        num_evaluated+=selection_pressure;
//
//        savefitness(num_evaluated,best_fitness,file_name_fitness);
//        //savefitnesshistory(i_generation,best_fitness_arr);
//    }
//    //std::cout<<best_position_history<<std::endl;
//
//
////    best_i = orderFitness(best_fitness_arr)(0);
//
//
//    std::cout<<"springCoefa<<"<<springCoefa_arr.row(best_i).transpose().format(print_formating)<<std::endl;
//    std::cout<<"springCoefb<<"<<springCoefb_arr.row(best_i).transpose().format(print_formating)<<std::endl<<std::endl;
//
//
//
//
//
//

    //////////////////////////////////////////////////////////////// playback ////////////////////////////////////

    InitialHeight = 0.5;

    int num_frames = 800;
    int skip_frames = 32;
    ArrayX3dRowMajor best_position_history = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);
    ArrayXdRowMajor springCoefa(num_cubes);
    ArrayXdRowMajor springCoefb(num_cubes);


//    // 1024
//    best=15.3847
//    springCoefa<<-0.0787365553522187,0.0397860520993408,0.0285811549930746,0.0148356401616873,0.0225603931014757,0.00942131270161891,0.0190325777693803,0.0233728654418946,0.0583954952556619,-0.0483040296883807,0.0242124729536801,0.00489416932481162,0.0665881243375549,0.0251908738910303,0.0741066694232201,0.0670986463979148,-0.0345643036600967,-0.0259200940029323,-0.0768512000128865,-0.066696290442113,0.0769517727367604,0.0440311868354987,0.00310712193810026,0.078899163230741,0.0788787391934348,0.0771249259249889,0.0117853437113109,0.0767958074048142,0.0437036912719271,0.0292698902400536,0.0590928833160164,0.0449671359942142,-0.0638690479769693,0.0677258177135295,0.0741325108270997,0.0139300205356702,-0.0556104941552554,0.078835353850776,0.00860470357830995,-0.0199256604630154,0.0330242166540873,0.0383049539282475,0.0671731600105638,0.0100578547688471,-0.0763951291127108,0.0609859394193561,0.0691503421495137,0.0247812721674833,0.0637547364010265,-0.0723889178002202,-0.0519388390946849,0.0624782460846371,-0.0431796303219998,0.00900752563448974,-0.0400732179545207,0.0376800872188434,0.021020399800878,0.0117542496590411,0.0336049248422409,0.0625928757258658,0.0278127577633072,0.0702474944993143,0.0435802303085012,-0.0415860856145556,0.0662081225558886,-0.056616993684609,0.0702688369854674,0.0669085326922944,-0.0577816398338329,0.00498168028038696,0.064275964950129,0.0082106890540437,-0.0675601910741814,-0.0297229902212149,0.0547141934895084,-0.0238513961172902,0.0504024790089589,-0.0269918353608771,0.0652930952540101,0.0246078541616946,0.0783054675153948,0.0217986215452878,-0.0539794790623614,-0.0749792838073239,0.0455257733005685,0.0198429829005026,0.0347960266260412,-0.0795648236570716,-0.060654577007822,-0.0186328556102854,0.079048270806227,0.0673910252187208,-0.0283853611854768,0.0768474262419641,4.22997400361602e-05,-0.0458803351064587,0.0660115074301192,-0.00968886334900225,0.047408684808348,-0.0717701324782195,0.035807761529371,-0.0174382546997807,0.0189957893067148,0.0482475704551896,0.0328387551535101,0.0190584340176443,-0.0556038256621006,0.00324123408796323,0.0127268029229173,0.0472632276808765,-0.052150911824848,0.0281599258042537,0.0174305656260953,-0.0261303908872094,0.035257932212356,-0.017043661147842,0.0726981328223454,-0.0570304115009636,0.0253911636689594,0.0197565364348446,0.016967082249871,0.0102833213978636,-0.0188242315961161,0.0559513717032743,-0.0549317130329701,0.0612180681439201,0.0193526401982814,-0.0689202056773567,0.0065515928256988,0.00290657720210824,-0.0606903381555762,-0.0726630338247228,-0.0690636828490829,0.0193897279593573,0.0555845366304668,0.0437750723044272,0.0468450873002075,0.0394696528074704,-0.0329836936821154,0.0769758056090101,-0.0703300195887359,-0.00513460550696338,0.0571396772193301,0.0142195497380376,0.0164515214097471,0.0308472113656125,0.0694353171229428,0.0780693864068339,0.0725435366842717,0.00753803458415812,-0.0130345523790617,0.0691651284023265,-0.0578139640287561,-0.0610438296855631,-0.0769839821741842,-0.00349377188994259,-0.0491218169634798,0.0680368567239194,0.0104728014513025,0.0277090399468825,0.0643284656034449,0.0148744277499296,0.0319480187967178,0.0469556046268861,0.0225398461253102,0.0540060857562377,-0.0767194037869197,-0.0488786590699473,0.0358852142383299,0.0235482812711195,0.00388935815724048,-0.0344991294827774,0.0072645609730069,0.0182211471806379,-0.0330087505062151,0.066253237997802,-0.007020162309995,0.0507489917067276,0.0579286789791326,-0.0554916582887488,-0.0443672990633023,-0.0351058734744349,0.0638490229759595,0.062662030843652,-0.016149703159998,0.0502813330899371,0.054324965017999,0.0552459516260987,0.0749330744627788,0.075759350480195,0.058129147798243,0.0282569440958355,-0.0268196486434059,-0.022496843981788,0.0128301607010894,0.00652780320799341,0.00679419323244019,0.0513311464426714,-0.0737902242661408,0.00482908282653852,-0.028651634505322,-0.00741968561309375,-0.0113375247695192,-0.0309082336495203,-0.0691985385069617,0.0356537247242656,-0.0510489410772216,0.00378129925754912,-0.0162516089231947,-0.0731202621725948,0.0282896409688004,0.019381092013503,0.00626741023776652,0.0410187031411662,0.0771998289214446,0.0356241612674781,0.00990672205602387,-0.0794443978366649,0.0199252965021531,-0.0576569127932456,0.00239666443429732,0.000663004027198092,-0.0266708465417246,0.0555770157908914,0.0329926150670633,0.032585694171761,-0.017895181075621,0.070577879464054,0.0374795822822829,-0.0116854053417618,-0.00815244164697474,0.00814443644515446,0.0608949091196502,0.0605100336580118,0.05723620287014,0.0716963706126885,0.0161637583077716,-0.0738127382815875,-0.00452233020426814,0.0799121494590827,-0.0669330003796765,0.0132565554729351,0.0192932413980799,-0.0151591360267061,0.0587106682500125,0.0164930703195245,-0.0595349748337339,0.0651585560967953,-0.0119821358900434,0.0313351379853371,0.0102307653847294,0.0735335526123895,-0.0783080907903184,0.00670429681703182,0.0459915443351453,0.029551836526744,0.01614561294026,-0.0519036367404757,0.0165703120532307,0.0635226350869781,0.0793850834321845,-0.0715821296682498,0.00138656943146689,0.0789834328713299,0.068927903989762,-0.0416776768684749,0.0791554930915526,0.00509166229753367,-0.0354904150755566,0.0487799919627843,0.0385185939230934,0.0429319663492875,-0.0717527819386464,-0.0557029469943153,0.0424174485180608,0.0275694849097959,0.0283013731032796,0.0628824737588328,0.065324034305884,0.00756486963832979,-0.0580133232371944,0.000926739197192052,0.040198942275384,0.0671719744676738,-0.0315740670038266,0.00459441494531404,0.0406146016064966,-0.0549752013268765,0.0632378235381273,0.0168862972014782,0.0632471188079785,0.0675372461255414,0.00629654168699845,-0.0716245620658736,0.0733686844228621,-0.0731565964935145,0.0207637777727245,-0.063657464228411,0.0332426104476874,0.0721893606780537,-0.059819667255422,0.0611593644201538,0.0371006349274426,0.008756611388529,-0.0662300293828501,0.0403580811924986,0.0173801001790549,-0.0183408645998411,0.0773245660948216,0.0381073966427275,-0.0354915834756063,0.0505520611814888,0.0779140615453543,0.0467182230794422,0.00691535268300928,-0.00752095301054463,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.072972533086243,-0.0563313566038065,0.00493563063820985,0.0168609884832338,-0.0755998318621888,-0.0764228115959199,0.0651290465822113,0.0677687595727708,-0.0255279210701249,0.0154125595537073,0.0310159908752032,0.0324546150026772,0.0672664034017367,-0.014407729028914,0.0427632872586899,-0.0203874216696189,-0.0329939846475581,0.0582008473473605,0.049530052081463,0.0420176874102258,-0.0526738212130376,0.0566280047414426,-0.0299665013467737,0.0202939150204388,-0.00878208395502627,0.0613053419223006,0.0651859232358668,0.0415810448984855,0.0684167652523223,0.0225791594861915,-0.0764551426453773,-0.00469409365425543,0.0657973458551789,-0.0340542982304722,-0.0545746671290019,-0.0363075187598856,-0.0704504885480043,0.0251474474115052,0.0239696648456015,0.00888312433328625,-0.0184900807302865,-0.026193907384851,0.0346668059912821,-0.0502069033543611,0.0727939190186036,-0.0652426554193919,-0.073405683857112,0.058832445019454,0.0433435391706107,-0.0276349627727805,-0.0627231468924895,0.0216381766375332,0.0484023522030692,0.0511634874649726,0.0491580558268187,-0.0303264040268615,0.0372968669781912,0.0201023868931081,0.0588998787006828,0.0257136321560077,-0.0356257363947228,0.0141284905145072,-0.0589804615727535,-0.0498283906140497,-0.0516095621751666,-0.0335551287017554,-0.00613590937393527,-0.0420600507231709,0.0715923187097499,-0.0621662446028396,-0.0562829418369955,-0.0100930846902044,-0.00836015198769055,0.00148987960139749,0.00289533462510226,0.0398096824344578,0.0162472241820056,0.00948965076799023,-0.0518638096618763,0.0600174140898688,0.0472304060902588,-0.0355273761393164,0.0202883375903072,0.0234392278378081,0.0422216261339447,-0.0518900550770993,0.0370514884947964,0.0149342752690121,0.0455472842443489,0.0482517623753528,0.0192839286938701,0.05911724067252,-0.00100677106576356,-0.0409349269237998,0.0659956001425141,0.0674249651969527,-0.0747426529204206,0.0431735729416369,0.0108531262403601,0.00681786573268896,4.67374548533606e-05,0.0628869602749529,-0.0176661403186927,0.0786480859660767,0.0524744645035985,-0.0418882346720846,0.00538064945739723,0.0144900607737209,0.0412227917798612,0.0745418645046977,0.0701532715513153,0.000597512655533077,-0.0725604402984308,0.0226924941421917,0.0577535415337205,0.0666303270992964,0.0287497000157599,-0.0541356280325612,0.0614123152109852,-0.0198186871129175,0.0215594750847957,0.0270599850426697,0.0423853336672181,-0.0427899707680521,0.0499074528313328,0.0793047787834378,0.0473748562227119,0.0524015135663986,-0.0379093903107147,0.0793937288558049;
//
//    springCoefb<<-2.43192959884493,6.08514030121006,4.04035678062355,4.0939778826623,0.885759964104128,-4.40632201143234,-3.84532036998163,-5.33669165449445,-0.201902543387136,-5.74421405700816,1.84616840882136,3.98087731087271,1.11369296464456,0.934977743915378,2.02614908115595,0.959073349058476,5.15053093739066,3.006446367785,-1.6732576622727,-1.23768408693475,1.27684178065925,6.04271727318241,3.77882575946868,-6.08085982355576,0.652931759876538,-5.48543869292053,1.81651248891278,5.96869830004477,1.31094790673134,3.23199520554859,0.163817214486977,5.16220350788641,-2.7309598098749,5.48964219698341,5.39277565112383,3.98484668574114,1.44712120543646,5.41086102056708,5.39994485533409,0.946468807561242,4.0832679471894,0.279088417230802,-4.91215907281787,4.6082018532571,0.92838679515286,-1.64664393517429,3.42329263641773,5.70322596816668,-2.65748894621183,1.29342632652648,0.568687125340165,-2.43705856028738,-1.92000280751018,-1.8801132043928,-0.517201690556528,1.68458566290145,4.00263656362903,2.50000791231667,0.452871199207279,-4.12575887624951,1.3802421085177,2.0898833298572,1.09336188408456,4.51445899470899,0.99940722926591,-2.21682713707008,1.5502626395456,5.9937160473333,-3.08915132235467,5.79470203052607,5.86985869733578,1.12879047483547,2.60477802497352,0.454742152576666,0.517470415219057,-5.5002405095043,-0.520598258586578,-5.43747913056799,-2.17154290659146,-4.77267751366817,0.0886738094516904,2.80624253330262,6.11985611039204,-5.33693277264827,-6.2580654340777,1.58892567487088,1.83628608125866,-5.35466976461089,2.3785831099318,-3.87135237064131,3.14285901158726,4.07157566778108,4.50171615921588,4.25724858523384,-4.71109470247867,-0.420077554397874,2.01939375860174,3.12235313706693,3.04526415542555,5.21342763624706,-3.62827223889698,-0.737268266315717,1.96095282746985,5.25969098607654,6.00065908626095,3.66813714062273,6.04263567657224,-0.803124378177302,1.71191742881433,0.288838383269688,0.707383314006201,1.04678416576555,-2.49262183233593,0.54405421854657,4.53851627650281,-2.46750206641364,2.16747008949913,-0.845862873766654,2.96766701262683,1.696825034699,1.26885943120368,1.67349275491913,-1.65457166366446,-0.215499085192086,-0.443091929617129,-0.0824811602914577,3.72336076311685,-4.70688337686707,4.86531505313461,1.28618349288566,-5.77664094062,-0.588400267973177,0.716840377031068,4.52819929759633,-1.61189448189663,0.434314257440346,4.15175172698769,2.17904391568207,5.91437508511471,0.300796283741831,2.01864875880472,0.338573193269244,5.38597896424107,5.14448955286199,3.82100555758673,4.50086791031566,3.6169715871486,3.34167081581541,5.25594889042253,3.6242449893384,2.869005561767,5.92319047857269,4.68127868871957,-3.14173898982101,-2.97035686401395,-0.230860614223697,-0.657414860877198,0.394065895226855,1.19275490160846,5.56181821803563,-2.91342920436725,1.82579686459939,-1.63544369410512,1.79662060059091,3.02465413262939,-1.00109808835509,4.72917326240044,4.64586125535466,2.46719717703016,5.95293752900116,-0.728069876461334,5.10586141337541,2.35273521470591,1.07370842850819,-1.93135481169837,4.37060868490869,-5.19487679371209,0.391702183144271,-1.21653276126522,-6.13711309553413,-3.61258757467371,-4.63071240534989,1.85101678249404,0.308870603002905,-1.4892661951709,5.49822140170814,0.837830493970494,5.16937915605619,5.61545234316673,3.98285231687726,4.89567875026041,-6.23992361757923,-5.34797394185981,-2.9636930323549,1.39715120233631,-2.12644323785466,3.40947531887144,1.80457371889052,4.09084134539069,-4.52002808382579,1.6512736377946,5.63214168372345,3.85473970142402,4.58615899606649,0.42266490637997,-4.35980031027435,0.0279554131341941,1.51097331266788,0.96973753156521,5.09460785186898,1.65704541713375,3.6403351568915,1.05442214408828,0.661927881806177,4.70902626508569,-1.02537075450348,0.192632458299389,-2.51630261939429,1.84969643455691,3.83759113021315,1.5746427306926,5.27170062187801,-4.70822915029589,2.50985398883279,2.6840672558054,-6.16522383252989,-5.89977445487354,3.05056480308843,2.21451166528256,4.47425209051715,5.5492140964325,3.75304063429704,3.82320856838893,3.12076859200485,2.05601442451186,-2.03731173108278,5.0441534817305,-4.19921536235395,5.75684678158511,-0.269294192555964,-5.38779271633664,5.43889699630295,-2.91214424151613,-4.92389726981756,3.36810651670152,-4.48630318228211,0.33391717567896,-4.53575025726799,1.0605535709224,4.47792557723166,-3.00405579583632,3.00532537434112,-0.591793793133674,2.71423009991725,5.14505028459612,-2.34202100170558,-1.31113837281043,5.52846103557426,-4.83914922729897,3.79489135958476,3.12888332563001,-5.57312033671814,1.83308250263569,3.19973161421444,3.83083346113838,-5.96668672335076,3.33493105350654,2.5918017370172,-3.88271688570471,3.38475512344542,4.01203987647953,2.89624094641009,-0.173391873853497,-0.589636697054897,-1.62803646600724,-1.48343160398213,0.848696330955331,4.98906590967172,4.47871139269655,1.19791773758934,-4.04150873496511,3.90776009537591,3.43991137256052,2.61626586932033,1.49626440893562,2.84893529190212,4.7158110969248,-3.53724639953568,4.50491096743946,1.49275078504112,-5.77939771606788,4.17997546079867,4.1568390815018,-6.27213179494185,-5.02828562505671,-3.57958680001224,3.34456559070511,3.46028085870612,0.799359763371303,1.47706369856494,2.56591515626386,3.63946459053563,5.86683490040684,5.73396412750752,4.27641316775727,4.53509332871182,5.61569700322845,2.9393886033297,-2.31855975768557,-5.79668449207325,-5.13523515189644,4.72865749192539,-1.49676813258366,4.1860763159865,3.37062530542857,1.21079109463131,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,5.12690620844648,3.61694724017608,3.05138930367418,-4.01275011252449,-2.72372981137355,-0.0180124155448347,4.12303341309968,-3.74873780023526,-0.887138111904024,-5.41990213502721,-5.9774984601979,0.552717896425435,0.921612258771997,3.17838284590814,1.86797359438216,-1.47430087116196,1.35138812294626,-3.90908982234295,0.382622065522161,3.58663702089682,-1.89477216410266,0.842265657489956,5.32264774649661,4.07279963344888,5.93452553491607,4.58997732369303,0.786062578409914,0.325200986907274,-5.9917002755706,0.588543628002016,-3.55057119553441,1.28408373976525,-1.25341481552174,-2.18636580089857,1.02549154472145,-2.97758908615075,-3.02601767758674,-0.79860588073586,-1.31513570862061,-5.10264541624672,-4.04774808632133,4.46998918241493,1.57982980509186,5.8374249517381,0.290319230405863,1.86919610233567,0.185816115052586,2.65571266392906,3.43118334743,1.54324804109866,-0.19298985557722,-5.47214491078457,2.00854924420448,2.55367598272938,2.67066286720342,-4.40899082454341,-4.04429312445419,5.6046481394415,-1.30713353515645,-3.75280819417312,5.43198984070525,0.316535709131635,3.81446074559214,-2.10461017481648,5.52229986841057,-1.44323290968642,1.20098593903276,-3.78690300917618,4.04134640957772,6.16903543626382,-0.65473740484626,-5.61395256815843,4.35583941282708,5.25665177966897,5.8310232749945,1.77763132207111,0.842662682004633,-0.266345809952913,1.76766768752563,2.41652541542875,4.29474875067796,5.02212468610824,-3.23113956818491,-4.46801644912111,3.88677049527466,-5.81971337358793,2.74959988627578,-3.37151591823711,2.76592180759093,-3.16280122963266,2.80490867323494,2.55004818904424,0.689441449147486,-0.0351284046447834,6.06388826542954,3.55993961157682,0.535486416593954,3.05170925752847,-1.00278238102184,5.46116454437677,5.81344950372511,-1.74189744386769,-0.332778030647898,-4.03931177144349,4.25502161714161,4.25874652459808,1.4159745097112,3.55494616590075,4.0525786915951,5.09901468640174,4.31990406550899,1.04140371424045,5.82363717143333,-2.31027463831456,-0.906303222026946,-0.757257994087879,-6.03794887918608,0.696102823104212,1.00660345861309,4.57369870271486,5.17365517040888,0.38104645055342,3.68186372776027,3.74661889802106,3.58529704079631,4.04925640652263,4.59815072532965,1.95431786572231,-5.33829054200515,1.63610028456171;

//
//    //337, best = 12.6906
//    springCoefa<<0.00645010438116739,-0.0452042573528384,-0.0528561693303549,0.0614960696462058,0.0587237226552985,0.0605971876937999,-0.0574423591687541,0.00162055362091426,0.0226982862980563,-0.0412959367601648,-0.0351734028547878,0.0597910659666132,0.0649722829391119,-0.0624649280228023,-0.0323109093086379,-0.0327793055552893,0.000505311377581812,0.049832757235817,0.0657287252627912,0.0202775308968417,-0.0502318214486501,0.072758553248252,0.0622989767172411,-0.041058367770658,-0.0641986652390094,0.0443994888497514,0.0753305620420649,0.0671439663260914,-0.0701050852938113,0.0263846494146126,0.0345452125158837,0.0746591436093017,0.0159853069186142,0.0112381603754518,-0.0314767897648163,-0.063013058892923,0.0406466706658931,0.00976108905382505,0.0503478702764715,-0.000596802048663054,0.0754078967847898,-0.0719830608330588,0.0279418782704392,0.0100578547688471,-0.0763951291127108,0.0609859394193561,-0.0178263663211029,0.067796497339288,0.0656016605335508,0.00866293032125707,-0.0486564285907226,0.00180728945965287,0.0578659887435378,-0.00522383693848916,0.0434141521730526,0.0756601472923812,-0.0697969719720059,-0.0534413669879741,-0.0400651946477896,-0.039408459718995,0.0538153286379236,0.0669229818819663,0.0045821657296753,-0.04839159755427,-0.0788843618328145,-0.0130147170708583,0.0347281700906941,-0.0612425526889239,0.0487222983542468,0.000600546082215817,-0.025359115500636,0.0797777333854594,-0.0358447358551643,0.0220513501755823,0.0328351535428479,0.0584371582659547,0.0483988474907348,0.0209474711590202,-0.0700425295485382,0.0739585553081513,-0.0503895985942285,-0.0386989581392607,-0.00423415530670162,-0.0102658141079665,0.0360772049222501,-0.0408200032081548,-0.0146056668900911,0.0462802329502442,-0.0142613701216231,0.0253291384621193,0.033267608466521,0.0798004311322237,0.0122521203440857,0.0342593110884816,-0.0485911664220463,0.0133677585112712,-0.0587554060568825,0.0661370036686477,0.0321252058223473,0.0699668922973643,0.0663888378377952,-0.0732339097527945,0.00828629412648507,0.046605987162549,0.0452604440437911,0.0712442745638264,-0.0272949249051953,0.0136592914600201,0.0184036247595431,-0.0173374543792277,0.00761784669366564,-0.00686234835854841,0.0709595197439747,-0.0766163086875418,0.0628718376079909,-0.0199592075962383,-0.0374363118956966,-0.0795413078365574,0.05294493957929,-0.0099366753967184,0.0379331902218671,0.0712711431003981,-0.0631979684453448,0.0125179405673716,0.0487215249094747,0.061298857977812,-0.0564282343426851,-0.00577347498655946,0.0454319094705544,-0.0524168368300501,-0.0151337902690907,-0.0781356664738318,0.0781029409161317,-0.0316597777379955,-0.0783601456127875,-0.0649374878336384,0.0548888318868768,-0.0307623483011324,-0.0501697162399859,-0.0421290921243509,-0.0792839760329965,0.0325275791437942,0.0699586538178654,-0.0771570660905713,-0.0238203955925165,-0.0254792333885465,0.0645156831388771,-0.0514078767464533,-0.0250205411505981,-0.0210381690883255,0.0186554479313341,-0.0670873510032368,-0.0465552956455179,0.0100968095334698,0.0355733549511024,0.00925781380816261,0.0752201241237876,-0.0166394266191122,-0.0765156612529027,0.000831012241929301,0.0109437365508376,0.0411950555467283,0.00269534576809749,0.0776356627291175,0.0500097646264786,-0.0693146277914358,0.0494778787202564,0.0115796025523821,0.0535728519286834,0.0539394733188392,-0.00349679883732311,-0.0166353404226877,0.0345605409678819,0.0220555843143051,-0.0212501866096864,0.0380845878918126,-0.0622068019873494,0.0772291018656018,-0.0685790633357033,0.0639410707314606,0.0502643645484539,0.0430886261738318,0.0263082978018108,-0.0439429772477332,-0.0268145643672042,0.0416108464270881,0.0453148366349353,-0.0714154615958293,-0.05502858026653,0.0487991753820326,0.00941555064609998,0.0230978417016193,-0.0428502761399608,-0.0678891036603083,-0.0350381663977346,0.0138377789874752,0.0635057832208855,0.0346857251015891,-0.000683305375596186,0.0661701948015646,-0.0478875381349994,0.0758198958615865,0.0431847738143347,0.0338875305952716,0.0178754801013858,0.0689522270490515,0.0397418925986869,0.0356686781140364,0.0357151095362916,0.0416943377078019,-0.0156460134385368,0.0735725044844338,0.055469166005067,-0.0634479441416673,0.0568009820472454,0.0119277532826773,0.0461261558374977,-0.0022123968239,-0.0791129107769173,0.0383166061147659,0.0178734747031114,-0.0427796691482792,0.0451157937036435,-0.0126576257183485,-0.0277317163104805,-0.0281253143903452,-0.0128802714929358,0.0210761391842161,0.0220082237158248,0.0120233691912253,0.0561449861415406,-0.00809438217808231,0.0786453260475049,-0.0119126403759758,-0.0233481096957569,-0.0388320614205823,0.047301308086129,0.035777097547323,0.000220454819603102,0.0786234311753062,-0.0387848637806181,0.00075965311413615,0.0540997214007366,0.0606307568057646,0.0536171465244224,0.0243289248386067,-0.066775535487931,0.0726062303082964,0.0729330518838805,0.00645560590712699,0.0648017433214941,-0.0625841674127542,-0.0248217000555348,0.0299175369506318,0.029980878477984,0.0274465836339847,-0.0782077774397134,0.0718779354504673,0.000800799077165132,-0.021507566543998,0.00390130456718676,-0.0553322911147644,0.0716716522687372,0.0572229261373854,0.0728941562062131,-0.052950058492343,0.0437145691941094,-0.0710342002618286,0.06282703905498,-0.0360649760607933,-0.0724107691610282,0.0715347042877101,0.0446946770533429,0.0742001225027256,-0.01329639632874,0.0183118235777653,0.0185290472668265,-7.19318166710092e-05,-0.0261754899221358,0.0613510288956347,0.0430038731185648,-0.0413737466006418,0.0787668615573863,-0.0385490345016816,0.06854379034999,0.00352506835177778,0.0688975491323031,0.0703360129102766,0.0133933155017766,-0.0425797281985077,-0.0311715537082271,0.0793043083694318,-0.0179120193132721,-0.0607735025048133,0.0158890933827534,-0.0368487085946131,-0.0645253226834002,0.0727586323371441,-0.0334753501198605,-0.00660638045827224,0.0304920743920338,-0.0369775748797588,0.078491745021675,-0.0403574947083171,-0.0313003325235566,0.0616713589570483,-0.0215847386892814,-0.000304401554309015,0.0661880928213653,-0.0318696097805489,-0.0730029912260375,0.0666616739270565,-0.0699727649940051,0.0124550792028035,0.05758461288064,0.0567268158095872,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.0658178208702327,-0.0563313566038065,0.0500588696310571,0.00382756894749756,0.0518925429749734,0.0558042647272125,0.0707498069809516,-0.0618310180221829,0.018713348702865,-0.0787391865992635,0.0290553043545481,-0.032300965614757,0.060290464004637,-0.0729320703413952,0.0403569689744523,0.0704091845408125,-0.0749069848632938,0.00644400089667897,-0.00989521708800234,0.0712811080325772,-0.0765668570234286,-0.00289820831403985,0.0579427818851279,0.0246742094588996,0.00888519283797834,0.0708139018094795,0.0666861274031392,0.0407490727122636,0.067646947472727,-0.0343162461902556,0.0547388540090708,0.0588098871609242,0.0757598622868582,0.0152984704521012,0.00793581355732671,0.0608907974422401,-0.0271419077120451,-0.0165874863120669,0.0629400009582471,0.0264471505332958,-0.0679513868633431,0.0416363874474283,0.0738929074497386,-0.0587957865646089,-0.0326140969771026,0.0538741387461368,-0.0751046355790946,-0.0621736521935899,0.0624445047147779,-0.0164819216432431,0.0584281103212517,0.0647932145821144,0.07246049503042,0.0403427731806146,0.0597265232259997,-0.0272246888779219,0.0556399539744667,0.0464126505546331,-0.0664756162401641,0.0114652567782743,-0.0679035956356226,0.0682632378434126,-0.00972485613530728,-0.0721437334232702,0.00356170822100793,0.0782109574220194,0.0687470640189699,0.0711567967754381,0.0393813721441467,0.0516870649027112,0.00286695096775281,-0.00632791582789641,0.0041841929823126,0.015323839190916,0.0148762976820005,0.00604671576095394,0.066970552647007,0.0321053442408369,-0.00363004591485022,-0.0448411121800733,-0.0110662296093377,-0.0349048151610907,0.0475583955494493,0.0369585034050786,0.0076267851089959,0.0164150435926463,-0.0503057194362887,0.0302437751508522,0.0430163776702324,0.0711511828027485,0.0283171388266222,0.0130441321120803,-0.0471295148726225,-0.0636974783910892,-0.0660869613970104,-0.011920530019291,0.0495780079670148,-0.049571110219495,0.0136326548893157,0.0102290127101489,0.0409999842015095,0.0746587701116962,-0.0481900098026684,0.0243057803550296,-0.0231146314475288,-0.0344240333672725,0.0388808117740176,0.0163504530193053,-0.071885456960595,0.000757831372673543,0.00110938649676246,0.0241672417261485,0.0141816137957746,0.0226031510264627,0.0329475054670812,0.0382672214387096,0.00347060505462373,0.0484186518045229,-0.0534132137025768,0.0135126887670402,0.00535699721675226,-0.035495290139455,-0.0313867033791666,0.0747957561099328,-0.0678140540364264,0.0133154983078207,0.038102190922062,0.0497354563184713,-0.000275638140866361,-0.0282651541886223;
//    springCoefb<<-3.86450119769613,3.67396148701663,-1.13984006050458,2.03765259098959,2.64542699241669,2.95856247417384,5.80200878818607,2.40355156331383,-4.69143756620047,5.16828797133064,-5.37462081282968,6.23123557546584,-2.54123569356608,-0.853094411257804,6.18780783481239,-0.25840404330108,-1.43313828941817,3.5426954547073,3.0099828551495,3.1224356846546,2.45892734778446,-3.8763295338547,0.926145971748672,0.0116756982145462,2.95324370712336,-5.44658419795381,0.871890565935848,0.216047562429265,-2.90890836040524,1.55194944318959,3.2013099083262,0.443808215258561,1.50866681169399,3.96182480405718,-4.44289417283143,2.07789435944452,0.400524817803756,3.56156850061099,-5.81259077082683,-5.57503069118609,-1.07917516017452,2.16978918383382,4.43453714645039,4.6082018532571,0.92838679515286,-1.64664393517429,-0.667023388585808,-2.67701341922912,1.97550020128691,-0.931188792323655,4.7886848273943,4.8221907184018,0.422466080091317,-2.26725924159576,0.546693422196206,-0.0183688877131311,-1.37746925424591,-3.18839334253656,4.88263154649922,4.586231539247,5.84782613461427,0.626629874566872,0.470299461674814,-4.64607910140807,2.86791570537105,-5.18858848726607,4.55014480331692,2.9092112956244,-2.31733057297406,2.02596680108756,-3.71580874952019,3.29165047590883,-5.62185617605571,5.69597240246917,1.07067794109222,2.73814272982173,0.708732434591533,5.82309725519261,-2.13598771901082,5.23093968992955,-1.39127674298272,-3.63048809746819,3.76994520833135,-0.132346752726178,0.385437866787722,-1.96654657532412,6.13246956541236,5.29115381254182,1.12824528213933,4.73191590605991,6.22595551861379,-0.877647887452685,-0.924639425224893,-0.412395970647958,0.759458216990917,-4.33990892570552,0.682200742085968,-0.973582185543834,4.85248756991888,4.6480553691119,3.01167867561052,-5.14650637960131,6.11865829406152,1.67234227541926,2.65551944889275,2.29173604351676,4.21403339977311,-2.91893332236739,1.37190715137686,-4.20513952508938,-3.97117883243784,0.87583388661951,1.9510871980551,6.08195157589351,-5.53969807195834,5.11618064423015,-2.16778020528228,4.55860292152544,1.08789178226731,-0.527725255953519,-0.322841148056009,5.99484878800427,4.22137079401472,1.31345772080213,-0.485785197631894,0.48210202002668,-4.56082833002572,2.85335845942771,0.988860129289455,-2.24536976258571,3.84903450826306,-2.74962286407637,1.82795550024617,1.6256001749435,-3.44664833496869,-0.787624398580966,5.25927800430962,-3.37926496591605,3.72316823967611,2.58035980433033,2.54019413687263,4.6156947458934,5.43975856755787,2.60790289401073,-0.295271064611664,-4.88084775784971,3.49341678315288,4.02247423416613,5.96094036367573,2.78165760372606,-2.78843622763906,-0.645085984380277,1.78861006168353,6.14382046948618,0.256582355197096,1.51007647953178,1.26258164177793,-2.87941997592207,-1.91975026689218,1.60290827367695,1.15839546734389,1.68661444175206,5.13647060960058,2.30322758574664,3.2393792208596,-1.69980019792689,-1.46685573007976,2.530977220624,-1.68944269128416,-0.361290403018153,3.67082353949202,5.27181747078639,5.74493160075912,1.13392076050876,-0.556325836430818,3.85611990323198,1.8993823338411,1.73442295313914,-1.7780976905355,0.671389861171619,6.01458241540703,3.58570349446706,5.95180612610385,-2.4627082360909,3.44633875810156,-6.18204504949343,5.33055344344088,4.06038690834072,-2.77827981956384,-2.87238202345131,-0.619890023834005,2.30605028098078,-0.943097782080431,-1.76660461423343,-4.91671859633324,3.22169629640312,0.1680088172115,4.35319006237781,1.93313335724957,5.56927634149866,5.02831704853653,-0.679228309110092,1.00973872109804,5.45244241973347,-5.82849274860133,4.74144188213636,0.0550478928232125,2.35407478523977,-0.976663354292901,-4.13896229686305,5.43417379595952,5.47634982195368,0.0176421053672753,-4.31250758108131,-4.15203447617205,-1.84821157858276,-6.19067216984734,-1.61778971396461,1.76101648377523,5.46675390560844,3.82072013307404,-0.0835571657792664,-2.70419611612083,-4.35664292349462,0.48699159593178,3.49123019508241,-4.89456587399671,-3.07676073412252,2.04452018239291,4.37967688426079,5.81464396461367,0.31627829285278,-5.69736777529547,1.20160379954745,-4.15864061034301,3.57165551955702,2.84271682577515,2.64859942517223,5.17866501670691,5.25211324025017,-2.71949634303427,0.624163225027658,2.49824505047781,3.65708436270744,5.12630444678776,2.90130931449199,-1.05154782982085,1.56897320310851,5.07119542386594,0.555765350707364,-0.096620585134944,6.10906236158531,-0.0394390155789386,3.48236849874422,2.64517698704675,-5.83563261964716,0.690413493826636,-3.14214643590601,3.9723615268424,-2.50113191339377,-5.04565475749689,-3.09774938915601,1.7045850603693,0.636181815904972,5.42613998921852,6.222876169116,5.39500339261633,-2.52784178346676,2.58829038843656,4.29048320347157,-3.55891374906826,0.151453472198172,-1.36853877735244,5.51428064317356,0.833913712280088,-2.52541953056468,2.13240475181388,6.06555108245924,-2.58978894902704,-0.581028044578514,0.234138328036739,3.59677566583801,-0.755150888844877,-3.03782608294973,0.795958958730565,0.997268549245814,-2.59027350259689,-4.79681275329447,0.61445411167008,1.06370354363343,-1.01475946668824,5.37583775584291,4.24913935447742,-3.19961308722923,5.02841343461045,1.28503691865316,2.35071693187312,0.483530942805895,-1.30176782090766,-3.06574137409335,5.07251818715228,-3.80211525034611,4.25142331163446,4.36818761739734,-5.58602333247991,0.138433648143927,6.13340090585955,6.16294620370673,5.40702213048559,-5.61419983056808,2.72281371453249,-1.16344619192875,-3.59943430289109,5.42763223796726,5.20067858627806,4.5670087462405,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,-1.86624674250335,3.61694724017608,3.16654841230712,1.42595260414086,5.78247426784757,3.7389969328919,5.13597709712056,-3.76063306196901,1.0853808746539,-3.51597829789533,2.89331316825022,-0.332643038421456,3.42466613621231,-5.612904709214,2.05983730548252,-3.53340152319352,-3.22784769827422,0.577054459623963,-3.65364052533846,-4.1040107736403,4.63797248144646,5.35235838919404,1.01572823443095,5.44427954487404,-0.0628920087411143,4.81985223290786,-0.814363134953018,-1.13390855392051,5.87110518702955,4.07865340641918,-1.20756497996135,-5.67315378971611,3.60556657449736,-2.16231139627231,1.6799593601027,-6.03598319715729,-1.45293454920339,-1.23471358501071,-6.09307128039791,-4.58914961603963,-0.26229766407284,3.8522889637721,4.70418797143078,4.74331604025211,-5.13477907841529,1.71304583270003,2.79302841518877,4.86964934250228,-1.99078880107362,-1.4619341284941,3.83699946763499,2.90383938221128,5.60335877519937,-3.74520178685621,-0.0386473194376096,5.75302449854938,-1.33615078980585,5.43017474560937,-1.6640692612228,4.81472937403923,3.22564294617688,3.41155095881585,5.4247607901748,0.548024314822573,-5.03394563745646,0.821534944425819,0.795226317665279,2.10498929431094,4.15626489805758,0.985340237267368,1.49034060315219,-0.6754763105094,3.46419741967238,2.21797383795514,-2.21534547025729,0.823900802681173,-6.20416934440633,4.70791372758925,4.4274468758401,1.20449927133029,5.13511560223299,-1.2277971191069,0.518862301505801,5.6298397646089,1.96500419324311,0.491842784345617,-3.65519818196669,-4.21394764351673,3.31549764356342,5.20712532451287,5.06518178769495,-5.19550258909431,-0.580364002666663,-3.8491704484604,-4.20797030078416,2.56470039592437,1.28176379164859,-1.83821881252117,0.0627040354743924,5.22340504166907,-1.35097660994907,-0.0896522383052312,3.3624099848038,3.61650616453083,4.76752530938366,-5.48805359457544,3.44166457960994,-5.44256429224387,-4.53917744709877,-4.98690482860416,-6.2334293264194,0.405232031674293,0.093060658909597,1.54348373475234,6.20180004328616,6.24117558242967,1.04370269768405,-5.34840537980092,5.23027941649461,1.63820279940981,1.30411456968309,-0.795897492765567,3.33436945617133,3.18023088396248,4.76810862162557,0.430241927938388,5.93355923657904,4.91211572363006,5.46638395118735,-0.286921933798238;
//
//
//    //generation288
//    // best = 12.4546
//    springCoefa<<0.00645010438116739,-0.0452042573528384,-0.0528561693303549,0.0614960696462058,0.0587237226552985,0.0605971876937999,-0.0574423591687541,0.00162055362091426,0.0226982862980563,-0.0412959367601648,-0.0351734028547878,0.0597910659666132,0.0649722829391119,-0.0624649280228023,-0.0323109093086379,-0.0327793055552893,0.000505311377581812,0.049832757235817,0.0657287252627912,0.0202775308968417,-0.0502318214486501,0.072758553248252,0.0622989767172411,-0.041058367770658,-0.0641986652390094,0.0443994888497514,0.0753305620420649,0.0671439663260914,-0.0701050852938113,0.0263846494146126,0.0345452125158837,0.0746591436093017,0.0159853069186142,0.0112381603754518,-0.0314767897648163,-0.063013058892923,0.0406466706658931,0.00976108905382505,0.0503478702764715,-0.000596802048663054,0.0754078967847898,-0.0719830608330588,-0.0351002073823941,0.0100578547688471,-0.0763951291127108,0.0609859394193561,-0.0178263663211029,0.067796497339288,0.0656016605335508,0.00866293032125707,-0.0486564285907226,0.00180728945965287,0.0578659887435378,-0.00522383693848916,0.0434141521730526,0.0756601472923812,-0.0697969719720059,-0.0534413669879741,-0.0400651946477896,-0.039408459718995,0.0538153286379236,0.0669229818819663,0.0045821657296753,-0.04839159755427,-0.0788843618328145,-0.0130147170708583,0.0347281700906941,-0.0612425526889239,0.0487222983542468,0.000600546082215817,-0.025359115500636,0.0797777333854594,-0.0358447358551643,0.0034972490452663,0.0328351535428479,0.0584371582659547,0.0483988474907348,0.0209474711590202,-0.0700425295485382,0.0739585553081513,-0.0503895985942285,-0.0386989581392607,-0.00423415530670162,-0.0102658141079665,0.0360772049222501,-0.0408200032081548,-0.0146056668900911,0.0462802329502442,-0.0142613701216231,0.0253291384621193,0.033267608466521,0.0798004311322237,0.0122521203440857,0.0342593110884816,-0.0485911664220463,0.0133677585112712,-0.0587554060568825,0.0661370036686477,0.0321252058223473,0.0699668922973643,0.0663888378377952,-0.0732339097527945,0.0179371517932652,0.046605987162549,0.0452604440437911,0.0225797791511657,-0.0272949249051953,0.0136592914600201,0.0755267916052145,-0.0173374543792277,0.00761784669366564,-0.00686234835854841,0.0709595197439747,-0.0766163086875418,0.0628718376079909,-0.0199592075962383,-0.0374363118956966,-0.0795413078365574,0.05294493957929,-0.0099366753967184,0.0379331902218671,0.0712711431003981,-0.0631979684453448,0.00244068909550117,0.0487215249094747,0.061298857977812,-0.0564282343426851,-0.00577347498655946,0.0454319094705544,-0.0524168368300501,-0.0151337902690907,-0.0781356664738318,0.0781029409161317,0.0338453969261504,-0.0783601456127875,-0.0649374878336384,0.0548888318868768,-0.0307623483011324,-0.0501697162399859,-0.0421290921243509,-0.0792839760329965,0.0325275791437942,0.0699586538178654,-0.0771570660905713,-0.0238203955925165,-0.0254792333885465,0.0412558472433245,-0.0514078767464533,-0.0250205411505981,-0.0210381690883255,0.0186554479313341,-0.0670873510032368,-0.0465552956455179,0.0100968095334698,0.0355733549511024,0.00925781380816261,0.0752201241237876,-0.0166394266191122,-0.0765156612529027,0.000831012241929301,0.0109437365508376,0.0411950555467283,0.00269534576809749,0.00904667739246351,0.0500097646264786,-0.0693146277914358,0.0494778787202564,0.0115796025523821,0.0535728519286834,0.0539394733188392,-0.00349679883732311,-0.0166353404226877,0.0345605409678819,0.0220555843143051,-0.0212501866096864,0.0380845878918126,-0.0622068019873494,0.0772291018656018,-0.0685790633357033,0.0639410707314606,0.0502643645484539,0.0430886261738318,0.0263082978018108,-0.0439429772477332,-0.0268145643672042,0.0416108464270881,0.0453148366349353,-0.0714154615958293,-0.05502858026653,0.0487991753820326,0.00941555064609998,0.0230978417016193,-0.0428502761399608,-0.0678891036603083,-0.0350381663977346,0.0738404946000504,0.0635057832208855,0.0346857251015891,-0.000683305375596186,0.0794203804868315,-0.0478875381349994,0.0758198958615865,0.0665370003816378,0.0338875305952716,0.0178754801013858,0.0689522270490515,0.0397418925986869,0.0356686781140364,0.0357151095362916,0.0416943377078019,-0.0156460134385368,0.0779409802272265,0.055469166005067,-0.0634479441416673,0.0568009820472454,0.0119277532826773,0.0461261558374977,-0.0022123968239,-0.0791129107769173,0.0383166061147659,0.0178734747031114,-0.0427796691482792,0.0451157937036435,-0.0126576257183485,-0.0277317163104805,-0.0281253143903452,-0.0128802714929358,0.0210761391842161,0.0220082237158248,0.0120233691912253,0.0561449861415406,-0.00809438217808231,0.0786453260475049,-0.0119126403759758,-0.0233481096957569,-0.0388320614205823,0.047301308086129,0.035777097547323,0.000220454819603102,0.0786234311753062,-0.0387848637806181,0.00075965311413615,0.0540997214007366,0.0606307568057646,0.0536171465244224,0.0243289248386067,-0.066775535487931,0.037513221871031,0.0729330518838805,0.00645560590712699,0.0648017433214941,-0.0625841674127542,-0.0248217000555348,0.0299175369506318,0.0047582068688973,0.0274465836339847,-0.0782077774397134,0.0718779354504673,0.00371479947545656,-0.021507566543998,0.00390130456718676,-0.0553322911147644,0.0716716522687372,0.0572229261373854,0.0728941562062131,-0.052950058492343,0.0437145691941094,-0.0710342002618286,0.06282703905498,-0.0360649760607933,-0.0724107691610282,0.0715347042877101,0.0446946770533429,0.0742001225027256,-0.01329639632874,0.0183118235777653,0.0185290472668265,-7.19318166710092e-05,-0.0261754899221358,0.0613510288956347,0.0430038731185648,-0.0413737466006418,0.0787668615573863,-0.0385490345016816,0.06854379034999,0.00352506835177778,0.0688975491323031,0.0703360129102766,-0.00459699627226079,-0.0425797281985077,-0.0311715537082271,0.0793043083694318,-0.0179120193132721,-0.0607735025048133,0.0158890933827534,-0.0368487085946131,-0.0645253226834002,0.0727586323371441,-0.0334753501198605,-0.00660638045827224,0.0304920743920338,-0.0369775748797588,0.078491745021675,0.0263718062232503,-0.0313003325235566,0.0616713589570483,-0.0215847386892814,-0.000304401554309015,0.0661880928213653,-0.0318696097805489,-0.0730029912260375,0.00450031574266356,-0.0699727649940051,0.0124550792028035,0.05758461288064,0.0567268158095872,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.0658178208702327,-0.0563313566038065,0.0500588696310571,0.00382756894749756,0.0518925429749734,0.0558042647272125,0.0707498069809516,-0.0618310180221829,0.018713348702865,-0.0787391865992635,0.0290553043545481,-0.032300965614757,0.060290464004637,-0.0729320703413952,0.0403569689744523,0.0704091845408125,-0.0749069848632938,0.00644400089667897,-0.00989521708800234,0.0712811080325772,0.0746148275705074,-0.00289820831403985,0.0579427818851279,0.0246742094588996,0.00888519283797834,-0.0647028192806536,0.0666861274031392,0.0407490727122636,0.067646947472727,-0.0343162461902556,0.0547388540090708,0.0588098871609242,0.0757598622868582,0.0152984704521012,0.00793581355732671,0.0608907974422401,-0.0271419077120451,-0.0165874863120669,0.0629400009582471,0.0264471505332958,-0.0679513868633431,0.0416363874474283,0.0738929074497386,-0.0587957865646089,-0.0326140969771026,0.0538741387461368,-0.0751046355790946,-0.0621736521935899,0.0624445047147779,-0.0164819216432431,0.0584281103212517,0.0647932145821144,0.07246049503042,0.0403427731806146,0.0597265232259997,-0.0272246888779219,0.0556399539744667,0.0464126505546331,-0.0664756162401641,0.0114652567782743,-0.0679035956356226,0.0682632378434126,-0.00972485613530728,-0.0721437334232702,0.00356170822100793,0.0782109574220194,0.0687470640189699,0.0711567967754381,0.0393813721441467,0.0516870649027112,0.00286695096775281,-0.00632791582789641,0.0041841929823126,-0.040235586967429,0.0148762976820005,-0.0189104693843566,0.066970552647007,0.0321053442408369,-0.00363004591485022,-0.0448411121800733,-0.0110662296093377,-0.0349048151610907,0.0475583955494493,0.0369585034050786,0.0076267851089959,0.0164150435926463,-0.0503057194362887,0.0302437751508522,0.0430163776702324,0.0711511828027485,0.0283171388266222,0.0130441321120803,-0.0471295148726225,-0.0636974783910892,-0.0660869613970104,-0.011920530019291,0.0495780079670148,-0.049571110219495,0.0136326548893157,0.0102290127101489,0.0409999842015095,0.0746587701116962,-0.0481900098026684,0.0243057803550296,-0.0231146314475288,-0.0344240333672725,0.0388808117740176,0.0163504530193053,-0.071885456960595,0.000757831372673543,0.00110938649676246,0.0241672417261485,0.0141816137957746,0.0226031510264627,0.0329475054670812,0.0382672214387096,0.00347060505462373,0.0484186518045229,-0.0534132137025768,0.0135126887670402,0.00535699721675226,-0.035495290139455,-0.0313867033791666,0.0747957561099328,-0.0678140540364264,0.0133154983078207,0.038102190922062,0.0497354563184713,-0.000275638140866361,-0.0282651541886223;
//    springCoefb<<-3.86450119769613,3.67396148701663,-1.13984006050458,2.03765259098959,2.64542699241669,2.95856247417384,5.80200878818607,2.40355156331383,-4.69143756620047,5.16828797133064,-5.37462081282968,6.23123557546584,-2.54123569356608,-0.853094411257804,6.18780783481239,-0.25840404330108,-1.43313828941817,3.5426954547073,3.0099828551495,3.1224356846546,2.45892734778446,-3.8763295338547,0.926145971748672,0.0116756982145462,2.95324370712336,-5.44658419795381,0.871890565935848,0.216047562429265,-2.90890836040524,1.55194944318959,3.2013099083262,0.443808215258561,1.50866681169399,3.96182480405718,-4.44289417283143,2.07789435944452,0.400524817803756,3.56156850061099,-5.81259077082683,-5.57503069118609,-1.07917516017452,2.16978918383382,-4.66402037524789,4.6082018532571,0.92838679515286,-1.64664393517429,-0.667023388585808,-2.67701341922912,1.97550020128691,-0.931188792323655,4.7886848273943,4.8221907184018,0.422466080091317,-2.26725924159576,0.546693422196206,-0.0183688877131311,-1.37746925424591,-3.18839334253656,4.88263154649922,4.586231539247,5.84782613461427,0.626629874566872,0.470299461674814,-4.64607910140807,2.86791570537105,-5.18858848726607,4.55014480331692,2.9092112956244,-2.31733057297406,2.02596680108756,-3.71580874952019,3.29165047590883,-5.62185617605571,5.37772983878635,1.07067794109222,2.73814272982173,0.708732434591533,5.82309725519261,-2.13598771901082,5.23093968992955,-1.39127674298272,-3.63048809746819,3.76994520833135,-0.132346752726178,0.385437866787722,-1.96654657532412,6.13246956541236,5.29115381254182,1.12824528213933,4.73191590605991,6.22595551861379,-0.877647887452685,-0.924639425224893,-0.412395970647958,0.759458216990917,-4.33990892570552,0.682200742085968,-0.973582185543834,4.85248756991888,4.6480553691119,3.01167867561052,-5.14650637960131,5.25114733084475,1.67234227541926,2.65551944889275,-3.55598661973871,4.21403339977311,-2.91893332236739,1.66933174922785,-4.20513952508938,-3.97117883243784,0.87583388661951,1.9510871980551,6.08195157589351,-5.53969807195834,5.11618064423015,-2.16778020528228,4.55860292152544,1.08789178226731,-0.527725255953519,-0.322841148056009,5.99484878800427,4.22137079401472,3.86002698571005,-0.485785197631894,0.48210202002668,-4.56082833002572,2.85335845942771,0.988860129289455,-2.24536976258571,3.84903450826306,-2.74962286407637,1.82795550024617,4.92431787411572,-3.44664833496869,-0.787624398580966,5.25927800430962,-3.37926496591605,3.72316823967611,2.58035980433033,2.54019413687263,4.6156947458934,5.43975856755787,2.60790289401073,-0.295271064611664,-4.88084775784971,3.20907687217003,4.02247423416613,5.96094036367573,2.78165760372606,-2.78843622763906,-0.645085984380277,1.78861006168353,6.14382046948618,0.256582355197096,1.51007647953178,1.26258164177793,-2.87941997592207,-1.91975026689218,1.60290827367695,1.15839546734389,1.68661444175206,5.13647060960058,-3.29683423826162,3.2393792208596,-1.69980019792689,-1.46685573007976,2.530977220624,-1.68944269128416,-0.361290403018153,3.67082353949202,5.27181747078639,5.74493160075912,1.13392076050876,-0.556325836430818,3.85611990323198,1.8993823338411,1.73442295313914,-1.7780976905355,0.671389861171619,6.01458241540703,3.58570349446706,5.95180612610385,-2.4627082360909,3.44633875810156,-6.18204504949343,5.33055344344088,4.06038690834072,-2.77827981956384,-2.87238202345131,-0.619890023834005,2.30605028098078,-0.943097782080431,-1.76660461423343,-4.91671859633324,-3.67139857161773,0.1680088172115,4.35319006237781,1.93313335724957,0.0219546126774578,5.02831704853653,-0.679228309110092,2.52791093267738,5.45244241973347,-5.82849274860133,4.74144188213636,0.0550478928232125,2.35407478523977,-0.976663354292901,-4.13896229686305,5.43417379595952,5.50214937321627,0.0176421053672753,-4.31250758108131,-4.15203447617205,-1.84821157858276,-6.19067216984734,-1.61778971396461,1.76101648377523,5.46675390560844,3.82072013307404,-0.0835571657792664,-2.70419611612083,-4.35664292349462,0.48699159593178,3.49123019508241,-4.89456587399671,-3.07676073412252,2.04452018239291,4.37967688426079,5.81464396461367,0.31627829285278,-5.69736777529547,1.20160379954745,-4.15864061034301,3.57165551955702,2.84271682577515,2.64859942517223,5.17866501670691,5.25211324025017,-2.71949634303427,0.624163225027658,2.49824505047781,3.65708436270744,5.12630444678776,2.90130931449199,-1.05154782982085,0.384388370846745,5.07119542386594,0.555765350707364,-0.096620585134944,6.10906236158531,-0.0394390155789386,3.48236849874422,-4.53076576190931,-5.83563261964716,0.690413493826636,-3.14214643590601,2.9449940260152,-2.50113191339377,-5.04565475749689,-3.09774938915601,1.7045850603693,0.636181815904972,5.42613998921852,6.222876169116,5.39500339261633,-2.52784178346676,2.58829038843656,4.29048320347157,-3.55891374906826,0.151453472198172,-1.36853877735244,5.51428064317356,0.833913712280088,-2.52541953056468,2.13240475181388,6.06555108245924,-2.58978894902704,-0.581028044578514,0.234138328036739,3.59677566583801,-0.755150888844877,-3.03782608294973,0.795958958730565,0.997268549245814,-2.59027350259689,-4.79681275329447,4.1383073133398,1.06370354363343,-1.01475946668824,5.37583775584291,4.24913935447742,-3.19961308722923,5.02841343461045,1.28503691865316,2.35071693187312,0.483530942805895,-1.30176782090766,-3.06574137409335,5.07251818715228,-3.80211525034611,4.25142331163446,4.82321190907592,-5.58602333247991,0.138433648143927,6.13340090585955,6.16294620370673,5.40702213048559,-5.61419983056808,2.72281371453249,0.402073260036223,-3.59943430289109,5.42763223796726,5.20067858627806,4.5670087462405,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,-1.86624674250335,3.61694724017608,3.16654841230712,1.42595260414086,5.78247426784757,3.7389969328919,5.13597709712056,-3.76063306196901,1.0853808746539,-3.51597829789533,2.89331316825022,-0.332643038421456,3.42466613621231,-5.612904709214,2.05983730548252,-3.53340152319352,-3.22784769827422,0.577054459623963,-3.65364052533846,-4.1040107736403,1.61509414119198,5.35235838919404,1.01572823443095,5.44427954487404,-0.0628920087411143,-3.87413420294964,-0.814363134953018,-1.13390855392051,5.87110518702955,4.07865340641918,-1.20756497996135,-5.67315378971611,3.60556657449736,-2.16231139627231,1.6799593601027,-6.03598319715729,-1.45293454920339,-1.23471358501071,-6.09307128039791,-4.58914961603963,-0.26229766407284,3.8522889637721,4.70418797143078,4.74331604025211,-5.13477907841529,1.71304583270003,2.79302841518877,4.86964934250228,-1.99078880107362,-1.4619341284941,3.83699946763499,2.90383938221128,5.60335877519937,-3.74520178685621,-0.0386473194376096,5.75302449854938,-1.33615078980585,5.43017474560937,-1.6640692612228,4.81472937403923,3.22564294617688,3.41155095881585,5.4247607901748,0.548024314822573,-5.03394563745646,0.821534944425819,0.795226317665279,2.10498929431094,4.15626489805758,0.985340237267368,1.49034060315219,-0.6754763105094,3.46419741967238,5.63757408362887,-2.21534547025729,-1.06473449124515,-6.20416934440633,4.70791372758925,4.4274468758401,1.20449927133029,5.13511560223299,-1.2277971191069,0.518862301505801,5.6298397646089,1.96500419324311,0.491842784345617,-3.65519818196669,-4.21394764351673,3.31549764356342,5.20712532451287,5.06518178769495,-5.19550258909431,-0.580364002666663,-3.8491704484604,-4.20797030078416,2.56470039592437,1.28176379164859,-1.83821881252117,0.0627040354743924,5.22340504166907,-1.35097660994907,-0.0896522383052312,3.3624099848038,3.61650616453083,4.76752530938366,-5.48805359457544,3.44166457960994,-5.44256429224387,-4.53917744709877,-4.98690482860416,-6.2334293264194,0.405232031674293,0.093060658909597,1.54348373475234,6.20180004328616,6.24117558242967,1.04370269768405,-5.34840537980092,5.23027941649461,1.63820279940981,1.30411456968309,-0.795897492765567,3.33436945617133,3.18023088396248,4.76810862162557,0.430241927938388,5.93355923657904,4.91211572363006,5.46638395118735,-0.286921933798238;
//
////
////
////    generation180
////    best = 11.7085
//    springCoefa<<0.00645010438116739,-0.0452042573528384,-0.0528561693303549,0.0614960696462058,0.0215831993248236,0.0605971876937999,-0.0574423591687541,0.00162055362091426,0.0226982862980563,-0.0412959367601648,-0.0351734028547878,0.0597910659666132,0.0649722829391119,-0.0624649280228023,-0.0323109093086379,-0.0327793055552893,0.000505311377581812,0.049832757235817,0.0657287252627912,0.0202775308968417,-0.0502318214486501,0.072758553248252,0.0622989767172411,-0.041058367770658,-0.0641986652390094,0.0443994888497514,0.0753305620420649,0.0671439663260914,-0.0701050852938113,0.0263846494146126,0.0345452125158837,0.0746591436093017,0.0159853069186142,0.0112381603754518,-0.0314767897648163,-0.063013058892923,0.0406466706658931,0.00976108905382505,0.0503478702764715,-0.000596802048663054,0.0754078967847898,-0.0719830608330588,-0.0351002073823941,0.0100578547688471,-0.0763951291127108,0.0609859394193561,-0.0178263663211029,0.067796497339288,0.0656016605335508,0.00866293032125707,-0.0486564285907226,0.00180728945965287,0.0578659887435378,-0.00522383693848916,0.0434141521730526,0.0756601472923812,-0.0697969719720059,-0.0534413669879741,-0.0400651946477896,-0.039408459718995,0.0538153286379236,0.0669229818819663,0.0273875378572324,-0.04839159755427,-0.0788843618328145,-0.0130147170708583,0.0347281700906941,-0.0612425526889239,0.0487222983542468,0.000600546082215817,-0.025359115500636,0.0797777333854594,-0.0358447358551643,0.0384943537965856,0.0328351535428479,0.0584371582659547,0.0483988474907348,0.0209474711590202,-0.0700425295485382,0.0739585553081513,-0.0503895985942285,-0.0386989581392607,-0.00423415530670162,-0.0102658141079665,0.0360772049222501,-0.0408200032081548,-0.0146056668900911,0.0462802329502442,-0.0142613701216231,0.0253291384621193,-0.0731282267687508,0.0798004311322237,0.0122521203440857,0.0342593110884816,-0.0485911664220463,0.0133677585112712,-0.0587554060568825,0.0661370036686477,0.0321252058223473,0.0699668922973643,0.0663888378377952,-0.0732339097527945,0.0179371517932652,0.046605987162549,0.0452604440437911,0.0225797791511657,-0.0272949249051953,0.0136592914600201,0.0184036247595431,-0.0173374543792277,0.00761784669366564,-0.00686234835854841,0.0502887447576552,-0.0766163086875418,0.0628718376079909,-0.0199592075962383,-0.0374363118956966,-0.0795413078365574,0.05294493957929,-0.0099366753967184,0.0379331902218671,0.0712711431003981,-0.0631979684453448,0.00244068909550117,0.0487215249094747,0.061298857977812,-0.0564282343426851,-0.00577347498655946,0.0454319094705544,-0.0524168368300501,-0.0151337902690907,-0.0781356664738318,0.0781029409161317,-0.0316597777379955,-0.0783601456127875,-0.0649374878336384,0.0548888318868768,-0.0307623483011324,-0.0501697162399859,-0.0421290921243509,-0.0792839760329965,0.0325275791437942,0.0699586538178654,-0.0771570660905713,-0.0238203955925165,-0.0254792333885465,0.0412558472433245,-0.0514078767464533,-0.0250205411505981,-0.0210381690883255,0.0186554479313341,-0.0670873510032368,-0.0465552956455179,0.0100968095334698,0.0355733549511024,0.00925781380816261,0.0353991028458807,-0.0166394266191122,-0.0765156612529027,0.000831012241929301,0.0109437365508376,-0.0116494515219934,0.00269534576809749,0.00904667739246351,0.0366907707400111,-0.0693146277914358,0.0494778787202564,0.0115796025523821,0.0535728519286834,0.0539394733188392,-0.00349679883732311,-0.0166353404226877,0.0345605409678819,0.0220555843143051,-0.0212501866096864,0.0380845878918126,-0.0622068019873494,-0.00957170416115397,-0.0685790633357033,0.0413441979857923,0.00975890782184848,0.0430886261738318,0.0263082978018108,-0.0439429772477332,-0.0268145643672042,0.0416108464270881,0.0453148366349353,-0.0714154615958293,-0.05502858026653,0.0487991753820326,0.00941555064609998,0.0230978417016193,-0.0428502761399608,-0.0678891036603083,-0.0350381663977346,0.0738404946000504,0.0635057832208855,0.0346857251015891,-0.000683305375596186,0.00317234072981977,-0.0478875381349994,0.0758198958615865,0.0665370003816378,0.0338875305952716,0.0178754801013858,0.0689522270490515,0.0397418925986869,0.0356686781140364,0.0357151095362916,0.0416943377078019,-0.0156460134385368,-0.0345259826418599,0.055469166005067,-0.0634479441416673,0.0568009820472454,0.0119277532826773,0.0461261558374977,-0.0022123968239,-0.0791129107769173,0.0383166061147659,0.0178734747031114,-0.0427796691482792,0.0451157937036435,-0.0126576257183485,-0.0277317163104805,-0.0281253143903452,-0.0128802714929358,0.0210761391842161,0.0220082237158248,0.0120233691912253,0.0561449861415406,-0.00809438217808231,0.0786453260475049,-0.0119126403759758,-0.0233481096957569,-0.0388320614205823,-0.00378926869658254,0.035777097547323,0.000220454819603102,0.0786234311753062,-0.0387848637806181,0.00075965311413615,0.0666108916637538,-0.0373385715285961,0.0536171465244224,0.0243289248386067,-0.066775535487931,0.0355126865746047,0.0729330518838805,-0.0136554026294758,0.00602298338760805,-0.0625841674127542,-0.0248217000555348,0.0299175369506318,0.0047582068688973,0.0274465836339847,-0.0782077774397134,0.0666490129653202,-0.031477277256305,-0.021507566543998,0.012563850706896,-0.0553322911147644,0.0503980512779197,0.0572229261373854,0.0122499764406565,-0.052950058492343,0.0437145691941094,-0.0710342002618286,0.06282703905498,-0.0360649760607933,-0.0724107691610282,0.0715347042877101,0.0446946770533429,0.0742001225027256,-0.01329639632874,0.0183118235777653,0.0185290472668265,-7.19318166710092e-05,-0.0261754899221358,0.0613510288956347,0.066272665628359,-0.0413737466006418,0.0787668615573863,-0.0385490345016816,0.06854379034999,0.00352506835177778,0.0688975491323031,0.0703360129102766,-0.00459699627226079,-0.0425797281985077,-0.0311715537082271,0.0793043083694318,-0.0179120193132721,-0.0607735025048133,0.00185093883511188,-0.0368487085946131,-0.0645253226834002,0.0727586323371441,-0.0334753501198605,-0.00660638045827224,0.0304920743920338,-0.0369775748797588,-0.0693323536167538,-0.0403574947083171,-0.0313003325235566,-0.0441825657171116,-0.0215847386892814,-0.000304401554309015,0.0661880928213653,-0.0318696097805489,-0.0730029912260375,0.0666616739270565,-0.0699727649940051,0.0124550792028035,0.05758461288064,0.0567268158095872,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.0658178208702327,-0.0563313566038065,0.00166829061348676,0.00382756894749756,0.0518925429749734,0.0558042647272125,0.0707498069809516,-0.0618310180221829,0.018713348702865,-0.0787391865992635,0.0290553043545481,-0.032300965614757,0.060290464004637,-0.0729320703413952,-0.0231125085535983,0.0704091845408125,-0.0749069848632938,0.0353027527571203,-0.00989521708800234,0.0712811080325772,-0.0765668570234286,-0.00289820831403985,0.0579427818851279,0.0246742094588996,0.00888519283797834,-0.0647028192806536,0.0666861274031392,0.0407490727122636,0.0358253028783134,-0.0343162461902556,0.0547388540090708,0.0588098871609242,0.0757598622868582,0.0152984704521012,0.00793581355732671,0.0702418318562052,-0.0271419077120451,0.00573235547540119,0.0629400009582471,0.0264471505332958,-0.0679513868633431,0.0416363874474283,0.0738929074497386,-0.0587957865646089,-0.0326140969771026,0.0272061396144359,-0.0751046355790946,-0.0621736521935899,0.0535530162127939,-0.0164819216432431,0.0584281103212517,0.0647932145821144,0.07246049503042,0.0135721897701859,0.0597265232259997,-0.0272246888779219,0.0556399539744667,0.0464126505546331,-0.0664756162401641,0.0114652567782743,-0.0679035956356226,0.0682632378434126,-0.00972485613530728,-0.0721437334232702,0.00356170822100793,0.0782109574220194,0.0687470640189699,0.0711567967754381,0.0393813721441467,0.0516870649027112,0.00286695096775281,-0.00632791582789641,0.0041841929823126,-0.040235586967429,0.0148762976820005,0.0746886325996753,0.066970552647007,0.0321053442408369,-0.00363004591485022,-0.0448411121800733,-0.0110662296093377,-0.0349048151610907,0.0475583955494493,0.0369585034050786,-0.0023723154153546,0.0164150435926463,-0.0503057194362887,0.0302437751508522,0.0430163776702324,0.0711511828027485,0.0283171388266222,0.00998347676863824,-0.0471295148726225,-0.0636974783910892,-0.0660869613970104,-0.011920530019291,0.0495780079670148,-0.049571110219495,0.0136326548893157,0.0102290127101489,0.0409999842015095,0.00688007718561835,-0.0481900098026684,0.0243057803550296,-0.0231146314475288,-0.0344240333672725,0.0112182075955059,0.0163504530193053,-0.071885456960595,0.000757831372673543,0.00110938649676246,0.0241672417261485,0.0141816137957746,0.0226031510264627,0.0329475054670812,0.0382672214387096,0.00347060505462373,0.0484186518045229,-0.0534132137025768,0.0135126887670402,0.00535699721675226,-0.035495290139455,-0.0313867033791666,0.0747957561099328,-0.0678140540364264,-0.0456399912972189,0.038102190922062,0.0497354563184713,-0.000275638140866361,-0.0282651541886223;
//    springCoefb<<-3.86450119769613,3.67396148701663,-1.13984006050458,2.03765259098959,0.158803204430339,2.95856247417384,5.80200878818607,2.40355156331383,-4.69143756620047,5.16828797133064,-5.37462081282968,6.23123557546584,-2.54123569356608,-0.853094411257804,6.18780783481239,-0.25840404330108,-1.43313828941817,3.5426954547073,3.0099828551495,3.1224356846546,2.45892734778446,-3.8763295338547,0.926145971748672,0.0116756982145462,2.95324370712336,-5.44658419795381,0.871890565935848,0.216047562429265,-2.90890836040524,1.55194944318959,3.2013099083262,0.443808215258561,1.50866681169399,3.96182480405718,-4.44289417283143,2.07789435944452,0.400524817803756,3.56156850061099,-5.81259077082683,-5.57503069118609,-1.07917516017452,2.16978918383382,-4.66402037524789,4.6082018532571,0.92838679515286,-1.64664393517429,-0.667023388585808,-2.67701341922912,1.97550020128691,-0.931188792323655,4.7886848273943,4.8221907184018,0.422466080091317,-2.26725924159576,0.546693422196206,-0.0183688877131311,-1.37746925424591,-3.18839334253656,4.88263154649922,4.586231539247,5.84782613461427,0.626629874566872,2.27658908341489,-4.64607910140807,2.86791570537105,-5.18858848726607,4.55014480331692,2.9092112956244,-2.31733057297406,2.02596680108756,-3.71580874952019,3.29165047590883,-5.62185617605571,1.51884062849406,1.07067794109222,2.73814272982173,0.708732434591533,5.82309725519261,-2.13598771901082,5.23093968992955,-1.39127674298272,-3.63048809746819,3.76994520833135,-0.132346752726178,0.385437866787722,-1.96654657532412,6.13246956541236,5.29115381254182,1.12824528213933,4.73191590605991,3.59420015178882,-0.877647887452685,-0.924639425224893,-0.412395970647958,0.759458216990917,-4.33990892570552,0.682200742085968,-0.973582185543834,4.85248756991888,4.6480553691119,3.01167867561052,-5.14650637960131,5.25114733084475,1.67234227541926,2.65551944889275,-3.55598661973871,4.21403339977311,-2.91893332236739,1.37190715137686,-4.20513952508938,-3.97117883243784,0.87583388661951,3.63094354852836,6.08195157589351,-5.53969807195834,5.11618064423015,-2.16778020528228,4.55860292152544,1.08789178226731,-0.527725255953519,-0.322841148056009,5.99484878800427,4.22137079401472,3.86002698571005,-0.485785197631894,0.48210202002668,-4.56082833002572,2.85335845942771,0.988860129289455,-2.24536976258571,3.84903450826306,-2.74962286407637,1.82795550024617,1.6256001749435,-3.44664833496869,-0.787624398580966,5.25927800430962,-3.37926496591605,3.72316823967611,2.58035980433033,2.54019413687263,4.6156947458934,5.43975856755787,2.60790289401073,-0.295271064611664,-4.88084775784971,3.20907687217003,4.02247423416613,5.96094036367573,2.78165760372606,-2.78843622763906,-0.645085984380277,1.78861006168353,6.14382046948618,0.256582355197096,1.51007647953178,-5.66913705561251,-2.87941997592207,-1.91975026689218,1.60290827367695,1.15839546734389,-4.35390095862913,5.13647060960058,-3.29683423826162,3.55488441631438,-1.69980019792689,-1.46685573007976,2.530977220624,-1.68944269128416,-0.361290403018153,3.67082353949202,5.27181747078639,5.74493160075912,1.13392076050876,-0.556325836430818,3.85611990323198,1.8993823338411,-2.94824845053947,-1.7780976905355,1.41490136140737,4.45337492739446,3.58570349446706,5.95180612610385,-2.4627082360909,3.44633875810156,-6.18204504949343,5.33055344344088,4.06038690834072,-2.77827981956384,-2.87238202345131,-0.619890023834005,2.30605028098078,-0.943097782080431,-1.76660461423343,-4.91671859633324,-3.67139857161773,0.1680088172115,4.35319006237781,1.93313335724957,3.53927866774267,5.02831704853653,-0.679228309110092,2.52791093267738,5.45244241973347,-5.82849274860133,4.74144188213636,0.0550478928232125,2.35407478523977,-0.976663354292901,-4.13896229686305,5.43417379595952,-2.80647363275011,0.0176421053672753,-4.31250758108131,-4.15203447617205,-1.84821157858276,-6.19067216984734,-1.61778971396461,1.76101648377523,5.46675390560844,3.82072013307404,-0.0835571657792664,-2.70419611612083,-4.35664292349462,0.48699159593178,3.49123019508241,-4.89456587399671,-3.07676073412252,2.04452018239291,4.37967688426079,5.81464396461367,0.31627829285278,-5.69736777529547,1.20160379954745,-4.15864061034301,3.57165551955702,-0.631696188006531,2.64859942517223,5.17866501670691,5.25211324025017,-2.71949634303427,0.624163225027658,2.79000919224182,0.965119666877801,5.12630444678776,2.90130931449199,-1.05154782982085,6.21881578738931,5.07119542386594,3.50243225602164,1.61091039875577,6.10906236158531,-0.0394390155789386,3.48236849874422,-4.53076576190931,-5.83563261964716,0.690413493826636,1.74536302452319,-2.62920815376968,-2.50113191339377,5.97106231912336,-3.09774938915601,4.09833157945901,0.636181815904972,1.40544802961955,6.222876169116,5.39500339261633,-2.52784178346676,2.58829038843656,4.29048320347157,-3.55891374906826,0.151453472198172,-1.36853877735244,5.51428064317356,0.833913712280088,-2.52541953056468,2.13240475181388,6.06555108245924,-2.58978894902704,-0.581028044578514,3.28479813262921,3.59677566583801,-0.755150888844877,-3.03782608294973,0.795958958730565,0.997268549245814,-2.59027350259689,-4.79681275329447,4.1383073133398,1.06370354363343,-1.01475946668824,5.37583775584291,4.24913935447742,-3.19961308722923,-5.36718477694946,1.28503691865316,2.35071693187312,0.483530942805895,-1.30176782090766,-3.06574137409335,5.07251818715228,-3.80211525034611,3.43397377093537,4.36818761739734,-5.58602333247991,-5.26630890551536,6.13340090585955,6.16294620370673,5.40702213048559,-5.61419983056808,2.72281371453249,-1.16344619192875,-3.59943430289109,5.42763223796726,5.20067858627806,4.5670087462405,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,-1.86624674250335,3.61694724017608,3.14465653737126,1.42595260414086,5.78247426784757,3.7389969328919,5.13597709712056,-3.76063306196901,1.0853808746539,-3.51597829789533,2.89331316825022,-0.332643038421456,3.42466613621231,-5.612904709214,4.11877141200666,-3.53340152319352,-3.22784769827422,3.96898711201454,-3.65364052533846,-4.1040107736403,4.63797248144646,5.35235838919404,1.01572823443095,5.44427954487404,-0.0628920087411143,-3.87413420294964,-0.814363134953018,-1.13390855392051,-0.132305036154921,4.07865340641918,-1.20756497996135,-5.67315378971611,3.60556657449736,-2.16231139627231,1.6799593601027,5.20266362573856,-1.45293454920339,1.79485245078019,-6.09307128039791,-4.58914961603963,-0.26229766407284,3.8522889637721,4.70418797143078,4.74331604025211,-5.13477907841529,-5.55855822218352,2.79302841518877,4.86964934250228,2.86156468724705,-1.4619341284941,3.83699946763499,2.90383938221128,5.60335877519937,2.76661834140825,-0.0386473194376096,5.75302449854938,-1.33615078980585,5.43017474560937,-1.6640692612228,4.81472937403923,3.22564294617688,3.41155095881585,5.4247607901748,0.548024314822573,-5.03394563745646,0.821534944425819,0.795226317665279,2.10498929431094,4.15626489805758,0.985340237267368,1.49034060315219,-0.6754763105094,3.46419741967238,5.63757408362887,-2.21534547025729,0.0254404862144594,-6.20416934440633,4.70791372758925,4.4274468758401,1.20449927133029,5.13511560223299,-1.2277971191069,0.518862301505801,5.6298397646089,1.14711262406966,0.491842784345617,-3.65519818196669,-4.21394764351673,3.31549764356342,5.20712532451287,5.06518178769495,5.95868234964683,-0.580364002666663,-3.8491704484604,-4.20797030078416,2.56470039592437,1.28176379164859,-1.83821881252117,0.0627040354743924,5.22340504166907,-1.35097660994907,2.25602074960245,3.3624099848038,3.61650616453083,4.76752530938366,-5.48805359457544,4.09222998817019,-5.44256429224387,-4.53917744709877,-4.98690482860416,-6.2334293264194,0.405232031674293,0.093060658909597,1.54348373475234,6.20180004328616,6.24117558242967,1.04370269768405,-5.34840537980092,5.23027941649461,1.63820279940981,1.30411456968309,-0.795897492765567,3.33436945617133,3.18023088396248,4.76810862162557,-5.11639605598791,5.93355923657904,4.91211572363006,5.46638395118735,-0.286921933798238;
////
////
////    generation76
////    best = 9.81247
//    springCoefa<<-0.024/3335295209352,-0.0764124855382426,-0.0528561693303549,0.0614960696462058,0.0215831993248236,0.0605971876937999,-0.0574423591687541,0.00162055362091426,0.0226982862980563,-0.0412959367601648,-0.0351734028547878,0.0597910659666132,0.0649722829391119,-0.0624649280228023,-0.0323109093086379,-0.0327793055552893,0.0502365616997826,-0.0382879246577099,0.0657287252627912,0.0264999055613298,-0.0502318214486501,0.072758553248252,-0.0289847769350674,0.0567459386567863,-0.0641986652390094,0.0443994888497514,-0.0501613280783227,0.0671439663260914,-0.0701050852938113,0.0655792865835034,0.0345452125158837,0.0746591436093017,0.0159853069186142,-0.0710501655894565,-0.0314767897648163,-0.063013058892923,0.0406466706658931,0.00976108905382505,0.0503478702764715,-0.000596802048663054,0.0754078967847898,-0.0719830608330588,-0.0351002073823941,0.0100578547688471,-0.0763951291127108,0.0609859394193561,-0.0178263663211029,0.067796497339288,-0.0544402921825835,0.00866293032125707,-0.0486564285907226,0.00180728945965287,0.0578659887435378,-0.00522383693848916,0.0434141521730526,0.0756601472923812,-0.0697969719720059,0.0581262234510336,-0.0400651946477896,-0.039408459718995,0.0538153286379236,0.0669229818819663,0.0273875378572324,-0.04839159755427,-0.0788843618328145,-0.0130147170708583,0.0347281700906941,-0.0612425526889239,0.0487222983542468,-0.0797481658308525,-0.025359115500636,0.0797777333854594,0.0334085333004206,0.0384943537965856,0.0328351535428479,-0.0578390268133204,0.0483988474907348,0.0209474711590202,-0.0700425295485382,0.0666180210035646,-0.0503895985942285,-0.0386989581392607,-0.00423415530670162,-0.0102658141079665,0.0360772049222501,-0.0408200032081548,-0.0146056668900911,0.0462802329502442,-0.0142613701216231,0.0253291384621193,-0.0731282267687508,0.0798004311322237,0.0122521203440857,0.0342593110884816,-0.0485911664220463,0.0133677585112712,-0.0587554060568825,0.0661370036686477,0.0321252058223473,0.0699668922973643,0.0663888378377952,-0.0732339097527945,0.0697446256083179,-0.049455898091875,0.0452604440437911,0.0225797791511657,-0.0272949249051953,0.0136592914600201,0.0184036247595431,-0.0173374543792277,0.00761784669366564,-0.00686234835854841,0.0502887447576552,-0.0766163086875418,0.0628718376079909,-0.0199592075962383,-0.0374363118956966,-0.0795413078365574,0.05294493957929,-0.0099366753967184,0.0379331902218671,0.0544828735173135,-0.0631979684453448,0.00244068909550117,0.0487215249094747,0.00851160366484505,-0.0564282343426851,-0.00577347498655946,0.0454319094705544,-0.0524168368300501,-0.0151337902690907,-0.0781356664738318,0.0781029409161317,-0.0316597777379955,-0.0783601456127875,-0.0649374878336384,0.0548888318868768,-0.0307623483011324,-0.0501697162399859,-0.0421290921243509,-0.0792839760329965,0.0457219773371341,0.0699586538178654,-0.0771570660905713,-0.0238203955925165,0.0660915006352432,0.0465186778300063,-0.0514078767464533,-0.0250205411505981,0.0326495117886485,0.0186554479313341,-0.0670873510032368,0.00897859056311242,0.0100968095334698,-0.0402111922764271,0.00925781380816261,0.0353991028458807,-0.0166394266191122,-0.0765156612529027,0.000831012241929301,0.0109437365508376,-0.0116494515219934,0.00269534576809749,0.00904667739246351,0.00565128818723977,-0.0693146277914358,0.0494778787202564,0.0115796025523821,0.0535728519286834,0.0539394733188392,-0.00349679883732311,-0.0166353404226877,0.0345605409678819,0.0220555843143051,-0.0212501866096864,0.0380845878918126,-0.0622068019873494,-0.00957170416115397,-0.0685790633357033,0.0286853085219326,0.00975890782184848,0.0430886261738318,0.0263082978018108,-0.0439429772477332,-0.0268145643672042,0.0416108464270881,0.0453148366349353,-0.0714154615958293,-0.05502858026653,0.0487991753820326,0.00941555064609998,0.0359151562843077,-0.0428502761399608,-0.0678891036603083,-0.0350381663977346,0.0738404946000504,0.0635057832208855,0.0346857251015891,-0.000683305375596186,0.00317234072981977,-0.0478875381349994,0.0758198958615865,0.0665370003816378,0.026991156016938,0.0178754801013858,0.0689522270490515,0.0302734010435051,0.0356686781140364,0.0357151095362916,0.0416943377078019,-0.0156460134385368,-0.0345259826418599,0.055469166005067,-0.0634479441416673,0.0568009820472454,0.0433671771953476,0.0461261558374977,-0.0022123968239,-0.0791129107769173,0.0383166061147659,-0.0512571138009695,-0.069024978200451,0.0451157937036435,0.0591721793814425,-0.0277317163104805,-0.0281253143903452,-0.0128802714929358,0.0210761391842161,-0.0232997891042846,0.0120233691912253,0.0561449861415406,-0.00809438217808231,0.0786453260475049,-0.0119126403759758,-0.0233481096957569,-0.0388320614205823,-0.00378926869658254,0.035777097547323,0.000220454819603102,0.0786234311753062,-0.0387848637806181,0.00075965311413615,0.0666108916637538,-0.0373385715285961,0.0536171465244224,0.0243289248386067,-0.066775535487931,0.0355126865746047,0.0729330518838805,-0.0136554026294758,0.0648017433214941,-0.0625841674127542,-0.0248217000555348,0.0299175369506318,0.0047582068688973,0.0274465836339847,-0.0782077774397134,0.0718779354504673,-0.031477277256305,-0.021507566543998,0.00390130456718676,-0.0553322911147644,0.0503980512779197,0.0572229261373854,0.0127550685092598,-0.052950058492343,0.0437145691941094,-0.0710342002618286,0.06282703905498,-0.0360649760607933,-0.0724107691610282,0.0715347042877101,0.0446946770533429,0.0742001225027256,-0.01329639632874,0.0297103880542045,0.0185290472668265,-7.19318166710092e-05,-0.0261754899221358,0.0613510288956347,0.066272665628359,-0.0413737466006418,0.0787668615573863,-0.0385490345016816,0.06854379034999,0.00352506835177778,0.0688975491323031,0.0703360129102766,-0.00459699627226079,-0.0425797281985077,-0.0311715537082271,0.0793043083694318,-0.0179120193132721,-0.0607735025048133,0.00185093883511188,-0.0368487085946131,-0.0645253226834002,-0.0242658496388541,-0.0334753501198605,-0.00660638045827224,0.0304920743920338,-0.0369775748797588,-0.0693323536167538,-0.0403574947083171,-0.0313003325235566,-0.0441825657171116,-0.0215847386892814,-0.000304401554309015,0.0661880928213653,-0.0318696097805489,-0.0730029912260375,0.0666616739270565,-0.0699727649940051,-0.0682165988479818,0.05758461288064,0.0534941478008118,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,0.0556833512750491,0.0658178208702327,-0.0563313566038065,0.0500588696310571,0.0759067577849639,0.0518925429749734,0.0558042647272125,0.0707498069809516,-0.0618310180221829,0.018713348702865,-0.0787391865992635,0.0290553043545481,-0.032300965614757,0.060290464004637,0.0286717615658385,-0.0231125085535983,0.0704091845408125,0.0121482197138999,0.0353027527571203,-0.00989521708800234,0.0712811080325772,-0.0765668570234286,-0.00289820831403985,0.0579427818851279,-0.0665396220174337,0.00888519283797834,-0.0647028192806536,0.0666861274031392,0.0407490727122636,0.0358253028783134,-0.0343162461902556,0.0547388540090708,0.0588098871609242,0.0757598622868582,0.0152984704521012,0.00793581355732671,0.0608907974422401,-0.0271419077120451,-0.0165874863120669,0.0629400009582471,0.0264471505332958,-0.0679513868633431,-0.0379834373099652,0.036897462139324,-0.0587957865646089,-0.0326140969771026,0.0272061396144359,-0.0751046355790946,-0.0621736521935899,0.0624445047147779,-0.0164819216432431,0.0584281103212517,0.0647932145821144,0.07246049503042,0.0403427731806146,0.0597265232259997,-0.0272246888779219,0.0556399539744667,0.0464126505546331,-0.0664756162401641,0.0114652567782743,-0.0679035956356226,0.0682632378434126,-0.00972485613530728,-0.0721437334232702,0.00356170822100793,0.0782109574220194,0.0687470640189699,0.0711567967754381,-0.0183765289645533,0.0516870649027112,0.00286695096775281,-0.00632791582789641,0.0041841929823126,-0.040235586967429,0.0148762976820005,-0.0189104693843566,0.066970552647007,0.0321053442408369,-0.00363004591485022,-0.0448411121800733,-0.0110662296093377,-0.0349048151610907,0.0475583955494493,0.0369585034050786,-0.0023723154153546,0.0164150435926463,-0.0503057194362887,0.0302437751508522,0.0430163776702324,0.0649657570500699,0.0283171388266222,0.00998347676863824,-0.0471295148726225,-0.0636974783910892,-0.0660869613970104,-0.011920530019291,0.0495780079670148,-0.049571110219495,0.0136326548893157,0.0102290127101489,0.0409999842015095,0.0746587701116962,-0.0481900098026684,0.0243057803550296,-0.0231146314475288,-0.0344240333672725,0.0260604282657576,0.0163504530193053,-0.071885456960595,0.000757831372673543,0.00110938649676246,0.0241672417261485,-0.0195672743486088,0.0580885130007548,0.0734863896637626,-0.0458597718765306,0.0205103219023488,0.036194733733402,-0.0408694547977622,-0.0618210253034816,0.0313831599575389,0.0334226692251035,0.0638206570520162,0.00935702149260651,-0.0678140540364264,-0.0456399912972189,0.038102190922062,0.0497354563184713,-0.000275638140866361,-0.0282651541886223;
//    springCoefb<<2.17522656936397,6.09687438750377,-1.13984006050458,2.03765259098959,0.158803204430339,2.95856247417384,5.80200878818607,2.40355156331383,-4.69143756620047,5.16828797133064,-5.37462081282968,6.23123557546584,-2.54123569356608,-0.853094411257804,6.18780783481239,-0.25840404330108,5.18062186268341,-5.54559962383232,3.0099828551495,-2.71864908526521,2.45892734778446,-3.8763295338547,-1.19987948729439,1.68762737072376,2.95324370712336,-5.44658419795381,-0.0489276776337942,0.216047562429265,-2.90890836040524,2.65757145479013,3.2013099083262,0.443808215258561,1.50866681169399,0.359087494197565,-4.44289417283143,2.07789435944452,0.400524817803756,3.56156850061099,-5.81259077082683,-5.57503069118609,-1.07917516017452,2.16978918383382,-4.66402037524789,4.6082018532571,0.92838679515286,-1.64664393517429,-0.667023388585808,-2.67701341922912,-1.76097794466199,-0.931188792323655,4.7886848273943,4.8221907184018,0.422466080091317,-2.26725924159576,0.546693422196206,-0.0183688877131311,-1.37746925424591,5.46667785152358,4.88263154649922,4.586231539247,5.84782613461427,0.626629874566872,2.27658908341489,-4.64607910140807,2.86791570537105,-5.18858848726607,4.55014480331692,2.9092112956244,-2.31733057297406,-2.29792433884564,-3.71580874952019,3.29165047590883,3.91130103551232,1.51884062849406,1.07067794109222,-5.74215950563337,0.708732434591533,5.82309725519261,-2.13598771901082,6.03994030014279,-1.39127674298272,-3.63048809746819,3.76994520833135,-0.132346752726178,0.385437866787722,-1.96654657532412,6.13246956541236,5.29115381254182,1.12824528213933,4.73191590605991,3.59420015178882,-0.877647887452685,-0.924639425224893,-0.412395970647958,0.759458216990917,-4.33990892570552,0.682200742085968,-0.973582185543834,4.85248756991888,4.6480553691119,3.01167867561052,-5.14650637960131,1.65652063916906,3.67300769955481,2.65551944889275,-3.55598661973871,4.21403339977311,-2.91893332236739,1.37190715137686,-4.20513952508938,-3.97117883243784,0.87583388661951,3.63094354852836,6.08195157589351,-5.53969807195834,5.11618064423015,-2.16778020528228,4.55860292152544,1.08789178226731,-0.527725255953519,-0.322841148056009,5.29013766380914,4.22137079401472,3.86002698571005,-0.485785197631894,-0.636499661885087,-4.56082833002572,2.85335845942771,0.988860129289455,-2.24536976258571,3.84903450826306,-2.74962286407637,1.82795550024617,1.6256001749435,-3.44664833496869,-0.787624398580966,5.25927800430962,-3.37926496591605,3.72316823967611,2.58035980433033,2.54019413687263,3.00822696142147,5.43975856755787,2.60790289401073,-0.295271064611664,3.23356327376998,-4.20844502939635,4.02247423416613,5.96094036367573,1.9325862357431,-2.78843622763906,-0.645085984380277,1.2082300465003,6.14382046948618,-4.60177684589636,1.51007647953178,-5.66913705561251,-2.87941997592207,-1.91975026689218,1.60290827367695,1.15839546734389,-4.35390095862913,5.13647060960058,-3.29683423826162,0.173097849969219,-1.69980019792689,-1.46685573007976,2.530977220624,-1.68944269128416,-0.361290403018153,3.67082353949202,5.27181747078639,5.74493160075912,1.13392076050876,-0.556325836430818,3.85611990323198,1.8993823338411,-2.94824845053947,-1.7780976905355,-3.20308619513192,4.45337492739446,3.58570349446706,5.95180612610385,-2.4627082360909,3.44633875810156,-6.18204504949343,5.33055344344088,4.06038690834072,-2.77827981956384,-2.87238202345131,-0.619890023834005,4.66330084778005,-0.943097782080431,-1.76660461423343,-4.91671859633324,-3.67139857161773,0.1680088172115,4.35319006237781,1.93313335724957,3.53927866774267,5.02831704853653,-0.679228309110092,2.52791093267738,-5.47350201305081,-5.82849274860133,4.74144188213636,3.92232059952412,2.35407478523977,-0.976663354292901,-4.13896229686305,5.43417379595952,-2.80647363275011,0.0176421053672753,-4.31250758108131,-4.15203447617205,1.31873549888256,-6.19067216984734,-1.61778971396461,1.76101648377523,5.46675390560844,5.2620595089243,0.736745442809552,-2.70419611612083,4.51027418052252,0.48699159593178,3.49123019508241,-4.89456587399671,-3.07676073412252,3.09163979863127,4.37967688426079,5.81464396461367,0.31627829285278,-5.69736777529547,1.20160379954745,-4.15864061034301,3.57165551955702,-0.631696188006531,2.64859942517223,5.17866501670691,5.25211324025017,-2.71949634303427,0.624163225027658,2.79000919224182,0.965119666877801,5.12630444678776,2.90130931449199,-1.05154782982085,6.21881578738931,5.07119542386594,3.50243225602164,-0.096620585134944,6.10906236158531,-0.0394390155789386,3.48236849874422,-4.53076576190931,-5.83563261964716,0.690413493826636,-3.14214643590601,-2.62920815376968,-2.50113191339377,-5.04565475749689,-3.09774938915601,4.09833157945901,0.636181815904972,4.38703961039144,6.222876169116,5.39500339261633,-2.52784178346676,2.58829038843656,4.29048320347157,-3.55891374906826,0.151453472198172,-1.36853877735244,5.51428064317356,0.833913712280088,2.0496970306988,2.13240475181388,6.06555108245924,-2.58978894902704,-0.581028044578514,3.28479813262921,3.59677566583801,-0.755150888844877,-3.03782608294973,0.795958958730565,0.997268549245814,-2.59027350259689,-4.79681275329447,4.1383073133398,1.06370354363343,-1.01475946668824,5.37583775584291,4.24913935447742,-3.19961308722923,-5.36718477694946,1.28503691865316,2.35071693187312,-1.98863273488915,-1.30176782090766,-3.06574137409335,5.07251818715228,-3.80211525034611,3.43397377093537,4.36818761739734,-5.58602333247991,-5.26630890551536,6.13340090585955,6.16294620370673,5.40702213048559,-5.61419983056808,2.72281371453249,-1.16344619192875,-3.59943430289109,0.867934802064849,5.20067858627806,1.1712486788054,6.04617452801232,-4.80287018786472,-2.17921914776464,4.80479190166537,-1.86624674250335,3.61694724017608,3.16654841230712,-3.50020186774831,5.78247426784757,3.7389969328919,5.13597709712056,-3.76063306196901,1.0853808746539,-3.51597829789533,2.89331316825022,-0.332643038421456,3.42466613621231,3.09708842241205,4.11877141200666,-3.53340152319352,4.42904933476599,3.96898711201454,-3.65364052533846,-4.1040107736403,4.63797248144646,5.35235838919404,1.01572823443095,-5.24464702144464,-0.0628920087411143,-3.87413420294964,-0.814363134953018,-1.13390855392051,-0.132305036154921,4.07865340641918,-1.20756497996135,-5.67315378971611,3.60556657449736,-2.16231139627231,1.6799593601027,-6.03598319715729,-1.45293454920339,-1.23471358501071,-6.09307128039791,-4.58914961603963,-0.26229766407284,3.08470434990278,-2.13595171952332,4.74331604025211,-5.13477907841529,-5.55855822218352,2.79302841518877,4.86964934250228,-1.99078880107362,-1.4619341284941,3.83699946763499,2.90383938221128,5.60335877519937,-3.74520178685621,-0.0386473194376096,5.75302449854938,-1.33615078980585,5.43017474560937,-1.6640692612228,4.81472937403923,3.22564294617688,3.41155095881585,5.4247607901748,0.548024314822573,-5.03394563745646,0.821534944425819,0.795226317665279,2.10498929431094,5.87000655941511,0.985340237267368,1.49034060315219,-0.6754763105094,3.46419741967238,5.63757408362887,-2.21534547025729,-1.06473449124515,-6.20416934440633,4.70791372758925,4.4274468758401,1.20449927133029,5.13511560223299,-1.2277971191069,0.518862301505801,5.6298397646089,1.14711262406966,0.491842784345617,-3.65519818196669,-4.21394764351673,3.31549764356342,4.30680508554124,5.06518178769495,5.95868234964683,-0.580364002666663,-3.8491704484604,-4.20797030078416,2.56470039592437,1.28176379164859,-1.83821881252117,0.0627040354743924,5.22340504166907,-1.35097660994907,-0.0896522383052312,3.3624099848038,3.61650616453083,4.76752530938366,-5.48805359457544,5.39996857145514,-5.44256429224387,-4.53917744709877,-4.98690482860416,-6.2334293264194,0.405232031674293,2.71988049946384,2.20475918445201,-6.27230685941966,4.09078162907568,-4.44713727302319,1.28425551560327,-5.30665293779401,-2.82562864873845,1.9588105827259,-3.95291330573697,-1.13347551262121,1.94854899669419,4.76810862162557,-5.11639605598791,5.93355923657904,4.91211572363006,5.46638395118735,-0.286921933798238;
////
//////    generation29
//////    best = 7.76742
//    springCoefa<<-0.0301900419360912,0.000902070936409076,-0.0509075279398391,0.0457725778807758,0.0642844398991598,-0.0332345113685515,-0.0378343773436893,0.0331732349624733,0.00263208266470212,0.06589283078252,-0.034435630195977,0.0638287369086541,-0.0656598015062789,-0.0526862487814791,0.0254299279793305,0.0423902228299483,0.0594847144090965,-0.0796003606168555,-0.0702811236634297,0.0384057364419129,-0.0149202112550476,-0.0343579301025523,-0.0365674083477666,-0.041058367770658,-0.0641986652390094,0.0443994888497514,-0.0501613280783227,0.0671439663260914,-0.0701050852938113,0.0655792865835034,0.0345452125158837,0.0746591436093017,0.0159853069186142,-0.0710501655894565,-0.0314767897648163,-0.063013058892923,0.0406466706658931,0.00976108905382505,0.0503478702764715,-0.000596802048663054,0.0754078967847898,0.0678884030077087,-0.0400750472955755,0.0330574202318943,0.0580057091163498,-0.0700955063058508,0.0681123176161723,0.067796497339288,-0.0544402921825835,0.00866293032125707,-0.0486564285907226,0.00180728945965287,-0.039876215513738,0.0709734343695331,0.0434141521730526,0.0756601472923812,-0.0697969719720059,-0.0534413669879741,-0.0400651946477896,-0.039408459718995,0.0140618012538468,0.0669229818819663,0.0273875378572324,-0.04839159755427,-0.0788843618328145,-0.0130147170708583,0.0347281700906941,-0.0612425526889239,0.0487222983542468,-0.0797481658308525,-0.025359115500636,0.0797777333854594,-0.0358447358551643,0.0384943537965856,0.0328351535428479,-0.0578390268133204,0.0483988474907348,0.0197937520573982,-0.0700425295485382,0.0739585553081513,-0.0503895985942285,-0.0386989581392607,-0.00423415530670162,-0.0102658141079665,0.0360772049222501,-0.0408200032081548,-0.0146056668900911,0.0462802329502442,-0.0142613701216231,0.0253291384621193,-0.0731282267687508,0.0798004311322237,0.0122521203440857,0.0342593110884816,-0.0485911664220463,0.0133677585112712,-0.0587554060568825,0.0661370036686477,0.0321252058223473,0.0699668922973643,0.0663888378377952,-0.0732339097527945,0.0697446256083179,-0.049455898091875,0.0452604440437911,0.0225797791511657,-0.0272949249051953,0.0136592914600201,-0.0364727497643199,-0.0173374543792277,0.00761784669366564,-0.00686234835854841,0.0239635874815116,-0.0766163086875418,0.0628718376079909,-0.0199592075962383,-0.0374363118956966,-0.0795413078365574,0.012443153156174,-0.0099366753967184,0.0379331902218671,0.0544828735173135,-0.0631979684453448,0.00244068909550117,0.0487215249094747,0.00851160366484505,-0.0564282343426851,-0.00577347498655946,0.0454319094705544,-0.0524168368300501,-0.0151337902690907,-0.0781356664738318,0.0781029409161317,-0.0316597777379955,-0.0783601456127875,-0.0649374878336384,0.0548888318868768,-0.0307623483011324,-0.0501697162399859,-0.0421290921243509,-0.0792839760329965,0.0457219773371341,0.0699586538178654,-0.0771570660905713,-0.0238203955925165,-0.0254792333885465,0.0465186778300063,-0.0514078767464533,-0.0250205411505981,-0.0210381690883255,0.0186554479313341,-0.0670873510032368,-0.0465552956455179,0.0100968095334698,-0.0402111922764271,0.00925781380816261,0.0353991028458807,-0.0166394266191122,-0.0765156612529027,0.000831012241929301,0.0109437365508376,-0.0116494515219934,0.00269534576809749,0.00904667739246351,0.0366907707400111,0.00433520015530997,0.0494778787202564,0.0264849184949346,0.0497930037089591,-0.00257326316208265,-0.00349679883732311,-0.0166353404226877,0.0345605409678819,0.0220555843143051,-0.0212501866096864,0.0380845878918126,-0.0622068019873494,-0.00957170416115397,-0.0685790633357033,0.0286853085219326,0.00975890782184848,-0.0136816326965027,0.0321980692968695,0.0113269646891053,-0.0160849516541161,0.0295740999791651,0.0209866211288546,-0.0110406640595946,0.072190450351774,0.0309552830974363,-0.0699120674235337,-0.0411784830508653,0.0335775252960518,0.0226939648868022,0.0652078464558385,0.052831776203975,0.0333793370953665,0.0346857251015891,-0.000683305375596186,0.00317234072981977,-0.0478875381349994,0.0758198958615865,0.0665370003816378,0.0666730029073884,0.0178754801013858,0.0689522270490515,0.0302734010435051,0.0356686781140364,0.0357151095362916,0.0416943377078019,-0.0156460134385368,-0.0345259826418599,0.055469166005067,-0.0634479441416673,0.0568009820472454,0.0119277532826773,0.0461261558374977,-0.0022123968239,-0.0791129107769173,0.0383166061147659,-0.0512571138009695,-0.069024978200451,0.0451157937036435,-0.0126576257183485,-0.0277317163104805,-0.0281253143903452,-0.0128802714929358,0.0210761391842161,-0.0232997891042846,0.0120233691912253,0.0561449861415406,-0.00809438217808231,0.0786453260475049,-0.0119126403759758,-0.0233481096957569,-0.0388320614205823,-0.00378926869658254,0.035777097547323,0.000220454819603102,0.0786234311753062,-0.0387848637806181,0.00075965311413615,0.0666108916637538,-0.0373385715285961,0.0536171465244224,0.0243289248386067,-0.066775535487931,0.0355126865746047,-0.0371780183711918,-0.0136554026294758,0.0648017433214941,-0.0625841674127542,-0.0248217000555348,0.0299175369506318,0.0047582068688973,0.0274465836339847,-0.0782077774397134,0.0718779354504673,-0.031477277256305,-0.021507566543998,0.00390130456718676,-0.0553322911147644,0.0503980512779197,0.0491050120609471,0.0127550685092598,-0.052950058492343,0.0437145691941094,-0.0710342002618286,0.06282703905498,-0.0360649760607933,-0.0724107691610282,-0.0559578248001439,0.0446946770533429,0.0742001225027256,-0.01329639632874,0.0183118235777653,0.0185290472668265,-7.19318166710092e-05,-0.0261754899221358,0.0613510288956347,0.066272665628359,-0.0413737466006418,0.0787668615573863,-0.0385490345016816,0.06854379034999,0.00352506835177778,0.0688975491323031,0.0703360129102766,-0.00459699627226079,-0.0425797281985077,-0.0311715537082271,0.0793043083694318,-0.0179120193132721,-0.0607735025048133,0.00185093883511188,-0.0368487085946131,-0.0645253226834002,-0.0242658496388541,-0.0334753501198605,-0.00660638045827224,0.0304920743920338,-0.0369775748797588,-0.0693323536167538,-0.0403574947083171,-0.0313003325235566,-0.0441825657171116,-0.0215847386892814,-0.000304401554309015,0.0661880928213653,-0.0318696097805489,-0.0730029912260375,0.0666616739270565,-0.0699727649940051,-0.0682165988479818,0.05758461288064,-0.0243255691716101,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.0658178208702327,-0.0563313566038065,0.0500588696310571,0.0759067577849639,0.0518925429749734,0.0558042647272125,0.0707498069809516,-0.0618310180221829,0.018713348702865,-0.0787391865992635,0.0290553043545481,-0.032300965614757,0.060290464004637,0.0346382099178798,0.0522100783755118,-0.0166486700329225,0.0353441417754368,-0.0446828684046319,-0.0207169189773113,-0.00943085972658864,-0.0104108149234256,0.00108280252715703,0.00565539803619282,-0.029452938413924,0.0176620010741344,-0.0140700556775881,0.0785136153542034,0.0336327575676296,-0.0745179910932286,-0.00535490181546421,-0.00731789369476861,0.0270669642589367,0.0757598622868582,0.0152984704521012,0.00793581355732671,0.0608907974422401,-0.0271419077120451,-0.0165874863120669,0.0629400009582471,0.0264471505332958,-0.0679513868633431,-0.0379834373099652,0.036897462139324,-0.0587957865646089,0.0062197053326567,0.0272061396144359,-0.0751046355790946,-0.0621736521935899,0.0624445047147779,-0.0164819216432431,0.0584281103212517,0.0647932145821144,-0.0138886550506058,0.0403427731806146,0.0597265232259997,-0.0272246888779219,0.0556399539744667,0.0464126505546331,-0.0664756162401641,0.0114652567782743,-0.0679035956356226,0.0682632378434126,-0.00972485613530728,-0.0721437334232702,0.0365938888474669,0.0782109574220194,0.0687470640189699,0.0564198005089629,0.0640296162590522,0.0516870649027112,0.00286695096775281,-0.00632791582789641,-0.0662963724817598,-0.040235586967429,0.0791346156890151,-0.0189104693843566,0.066970552647007,0.0321053442408369,-0.00108412157794652,0.0494150573617849,-0.0767102596148431,-0.0226560113312006,0.0568617127728005,-0.0105989146654489,-0.062313238150586,0.0365882359988001,0.0421763964566292,0.0733267158238807,0.011341673693665,0.055700780290971,-0.0743169126633168,-0.0547912517072592,0.0395568902602218,-0.0240700065829186,0.0369082464077083,0.0524183360358786,0.0508395326746812,0.00847840397018456,0.0392911277846404,-0.0567403040718009,-0.0192200914487336,0.0674314556212311,0.0146710545265447,0.00483120601849221,-0.0061529317899388,0.0509080148352813,-0.0575229979201793,0.0563938426861511,0.0680953267906305,0.077270823864858,-0.0658073097401333,-0.0690422282316919,-0.0472917724248449,0.0226031510264627,0.0329475054670812,-0.0338470014714855,0.00347060505462373,0.0484186518045229,-0.0534132137025768,-9.4295702918572e-05,0.00535699721675226,-0.035495290139455,-0.0313867033791666,0.0633975579512294,-0.010011606109334,0.0712317100871502,0.00559826966635802,0.0555343693996627,0.0639076539799141,-0.00211184341558806;
//    springCoefb<<3.27777721437627,5.21464289550784,-3.60343900010493,-6.00399287169975,2.83043068310564,3.27712620413846,-1.22042688528079,0.383728839102502,-0.497761487470872,-2.74989623664747,-1.12757742493719,-2.36381488194657,-1.53953235212409,-5.23480154197283,0.990716244607597,5.18175423150528,2.54157246315401,-2.587319672087,3.69428281735015,0.44708529214188,-2.5166620583604,2.33178977414935,-0.772758863138392,0.0116756982145462,2.95324370712336,-5.44658419795381,-0.0489276776337942,0.216047562429265,-2.90890836040524,2.65757145479013,3.2013099083262,0.443808215258561,1.50866681169399,0.359087494197565,-4.44289417283143,2.07789435944452,0.400524817803756,3.56156850061099,-5.81259077082683,-5.57503069118609,-1.07917516017452,-2.22874335364205,2.24618422173198,4.06221267103506,6.16288187627402,5.47307700609747,-1.53076589175128,-2.67701341922912,-1.76097794466199,-0.931188792323655,4.7886848273943,4.8221907184018,-5.02425520974346,-3.89318801426694,0.546693422196206,-0.0183688877131311,-1.37746925424591,-3.18839334253656,4.88263154649922,4.586231539247,4.27729203625966,0.626629874566872,2.27658908341489,-4.64607910140807,2.86791570537105,-5.18858848726607,4.55014480331692,2.9092112956244,-2.31733057297406,-2.29792433884564,-3.71580874952019,3.29165047590883,-5.62185617605571,1.51884062849406,1.07067794109222,-5.74215950563337,0.708732434591533,4.32035811844568,-2.13598771901082,5.23093968992955,-1.39127674298272,-3.63048809746819,3.76994520833135,-0.132346752726178,0.385437866787722,-1.96654657532412,6.13246956541236,5.29115381254182,1.12824528213933,4.73191590605991,3.59420015178882,-0.877647887452685,-0.924639425224893,-0.412395970647958,0.759458216990917,-4.33990892570552,0.682200742085968,-0.973582185543834,4.85248756991888,4.6480553691119,3.01167867561052,-5.14650637960131,1.65652063916906,3.67300769955481,2.65551944889275,-3.55598661973871,4.21403339977311,-2.91893332236739,-4.01607457039777,-4.20513952508938,-3.97117883243784,0.87583388661951,-1.55244242255757,6.08195157589351,-5.53969807195834,5.11618064423015,-2.16778020528228,4.55860292152544,0.706917427270735,-0.527725255953519,-0.322841148056009,5.29013766380914,4.22137079401472,3.86002698571005,-0.485785197631894,-0.636499661885087,-4.56082833002572,2.85335845942771,0.988860129289455,-2.24536976258571,3.84903450826306,-2.74962286407637,1.82795550024617,1.6256001749435,-3.44664833496869,-0.787624398580966,5.25927800430962,-3.37926496591605,3.72316823967611,2.58035980433033,2.54019413687263,3.00822696142147,5.43975856755787,2.60790289401073,-0.295271064611664,-4.88084775784971,-4.20844502939635,4.02247423416613,5.96094036367573,2.78165760372606,-2.78843622763906,-0.645085984380277,1.78861006168353,6.14382046948618,-4.60177684589636,1.51007647953178,-5.66913705561251,-2.87941997592207,-1.91975026689218,1.60290827367695,1.15839546734389,-4.35390095862913,5.13647060960058,-3.29683423826162,3.55488441631438,-4.59336292536811,-1.46685573007976,5.63568670282605,-6.17683572703639,-5.60805821969295,3.67082353949202,5.27181747078639,5.74493160075912,1.13392076050876,-0.556325836430818,3.85611990323198,1.8993823338411,-2.94824845053947,-1.7780976905355,-3.20308619513192,4.45337492739446,-0.0701528982201201,-3.46349617118915,4.93762436242973,2.36090361650041,4.40502061708564,-3.74894043779257,-2.67395714349368,5.37424087545578,0.596664017037242,5.258914164886,-1.49210189797368,2.7103409629707,-1.15045581314572,-0.46313940754243,2.58063186027519,3.43292918892738,4.35319006237781,1.93313335724957,3.53927866774267,5.02831704853653,-0.679228309110092,2.52791093267738,4.49006344344398,-5.82849274860133,4.74144188213636,3.92232059952412,2.35407478523977,-0.976663354292901,-4.13896229686305,5.43417379595952,-2.80647363275011,0.0176421053672753,-4.31250758108131,-4.15203447617205,-1.84821157858276,-6.19067216984734,-1.61778971396461,1.76101648377523,5.46675390560844,5.2620595089243,0.736745442809552,-2.70419611612083,-4.35664292349462,0.48699159593178,3.49123019508241,-4.89456587399671,-3.07676073412252,3.09163979863127,4.37967688426079,5.81464396461367,0.31627829285278,-5.69736777529547,1.20160379954745,-4.15864061034301,3.57165551955702,-0.631696188006531,2.64859942517223,5.17866501670691,5.25211324025017,-2.71949634303427,0.624163225027658,2.79000919224182,0.965119666877801,5.12630444678776,2.90130931449199,-1.05154782982085,6.21881578738931,3.56975240360761,3.50243225602164,-0.096620585134944,6.10906236158531,-0.0394390155789386,3.48236849874422,-4.53076576190931,-5.83563261964716,0.690413493826636,-3.14214643590601,-2.62920815376968,-2.50113191339377,-5.04565475749689,-3.09774938915601,4.09833157945901,5.87266575426221,4.38703961039144,6.222876169116,5.39500339261633,-2.52784178346676,2.58829038843656,4.29048320347157,-3.55891374906826,6.15197924540229,-1.36853877735244,5.51428064317356,0.833913712280088,-2.52541953056468,2.13240475181388,6.06555108245924,-2.58978894902704,-0.581028044578514,3.28479813262921,3.59677566583801,-0.755150888844877,-3.03782608294973,0.795958958730565,0.997268549245814,-2.59027350259689,-4.79681275329447,4.1383073133398,1.06370354363343,-1.01475946668824,5.37583775584291,4.24913935447742,-3.19961308722923,-5.36718477694946,1.28503691865316,2.35071693187312,-1.98863273488915,-1.30176782090766,-3.06574137409335,5.07251818715228,-3.80211525034611,3.43397377093537,4.36818761739734,-5.58602333247991,-5.26630890551536,6.13340090585955,6.16294620370673,5.40702213048559,-5.61419983056808,2.72281371453249,-1.16344619192875,-3.59943430289109,0.867934802064849,5.20067858627806,-6.1831742934951,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,-1.86624674250335,3.61694724017608,3.16654841230712,-3.50020186774831,5.78247426784757,3.7389969328919,5.13597709712056,-3.76063306196901,1.0853808746539,-3.51597829789533,2.89331316825022,-0.332643038421456,3.42466613621231,3.2140630252074,-0.81281266417827,3.55729249646634,-3.99732763832671,5.16759031741356,-6.01963792010209,6.08607945080758,2.40865155897879,-0.34897150947404,-0.689990183858242,4.3800432750135,-1.6897918664172,4.00166629047414,1.29512865989868,-2.61690149022051,1.5858511172895,-4.8513777300623,2.07190843816047,-1.17473630883154,3.60556657449736,-2.16231139627231,1.6799593601027,-6.03598319715729,-1.45293454920339,-1.23471358501071,-6.09307128039791,-4.58914961603963,-0.26229766407284,3.08470434990278,-2.13595171952332,4.74331604025211,3.45098937921172,-5.55855822218352,2.79302841518877,4.86964934250228,-1.99078880107362,-1.4619341284941,3.83699946763499,2.90383938221128,-4.02582465590906,-3.74520178685621,-0.0386473194376096,5.75302449854938,-1.33615078980585,5.43017474560937,-1.6640692612228,4.81472937403923,3.22564294617688,3.41155095881585,5.4247607901748,0.548024314822573,5.20920498269117,0.821534944425819,0.795226317665279,-0.203694986659853,3.96434777754175,0.985340237267368,1.49034060315219,-0.6754763105094,-2.21314061868153,5.63757408362887,5.2510218292457,-1.06473449124515,-6.20416934440633,4.70791372758925,-2.47827035459454,-1.91177294547994,-0.884251189414299,-4.92445609281122,-5.9599007614408,1.37310935467664,-2.38647267966743,0.284637119121588,0.842948653226024,2.56056173052672,5.14163171950606,5.46206459200322,-5.0282914416189,-3.72990417564026,-1.5051535554682,3.52138756614046,1.70753641288865,-1.9099825502925,5.36430706499181,5.2016018805556,2.29881549392972,1.43413112508871,-3.10014662030768,3.47833679606309,2.05592782561468,-0.736004088604331,-5.48585112463027,2.08016178295133,-6.07073525442908,5.49252503872056,-0.897913872670351,4.8649391267239,-0.386430778530009,-3.67537460148936,2.32622207187887,1.54348373475234,6.20180004328616,-0.911167348488062,1.04370269768405,-5.34840537980092,5.23027941649461,-2.35946564997576,1.30411456968309,-0.795897492765567,3.33436945617133,6.01092648885484,-5.68498525066764,-1.79991515301554,3.55102793025452,0.61982206896785,-0.00910750460729953,1.52609584770269;
////
//
////    generation14
////    best = 6.76067
//    springCoefa<<0.00645010438116739,-0.0452042573528384,-0.0116153831461516,0.0502180757234889,0.0239324127621634,0.00482616853193667,-0.00859522130740584,-0.00404057171383899,0.0110200334764179,0.0562440407910589,-0.000325039141031462,-0.0320709737167093,0.0220207565007828,-0.00482695210949841,0.062534640497777,-0.0483893871160175,0.0248435026150399,-0.0577308989957585,0.0712340298254201,-0.0360334036294526,-0.0149202112550476,-0.0343579301025523,-0.0074666636844476,-0.0139881281992365,-0.0278726067896339,-0.047397683024126,0.0399086870066396,0.0345520573642813,0.0383247365608461,-0.0669448237060312,0.0104844747905081,-0.0624772769876184,-0.00109922882220672,0.0788690916443565,0.0677407987358704,-0.0571668160600433,0.0036952601017874,-0.0208544226460412,0.0187926122261177,-0.0652847064963005,-0.0446103819294881,-0.0615324269149138,-0.017355680138504,0.0574103745712947,0.0136406209755878,-0.0348210397152328,-0.0709790125447228,-0.0415158764093723,-0.0544402921825835,0.00866293032125707,-0.0486564285907226,0.00180728945965287,-0.039876215513738,-0.00522383693848916,0.0434141521730526,0.0756601472923812,-0.0697969719720059,-0.0534413669879741,-0.0400651946477896,-0.039408459718995,0.0140618012538468,0.0669229818819663,0.0273875378572324,-0.04839159755427,-0.0788843618328145,-0.0130147170708583,0.0347281700906941,-0.0612425526889239,0.0487222983542468,-0.0797481658308525,-0.025359115500636,0.0797777333854594,-0.0358447358551643,0.0384943537965856,0.0328351535428479,-0.0578390268133204,0.0483988474907348,0.0209474711590202,-0.0700425295485382,0.0739585553081513,-0.0503895985942285,-0.0386989581392607,-0.00423415530670162,0.0408158631844995,0.0360772049222501,-0.0408200032081548,-0.0146056668900911,0.0462802329502442,-0.0142613701216231,0.0253291384621193,-0.0731282267687508,0.0798004311322237,0.0122521203440857,0.0342593110884816,-0.0485911664220463,0.0133677585112712,-0.0587554060568825,0.0661370036686477,0.0321252058223473,0.0699668922973643,0.0663888378377952,-0.0732339097527945,0.0697446256083179,-0.049455898091875,0.0452604440437911,0.0225797791511657,-0.0272949249051953,0.0136592914600201,-0.0364727497643199,-0.0173374543792277,0.00761784669366564,-0.00686234835854841,0.0239635874815116,-0.0766163086875418,0.0628718376079909,-0.0199592075962383,-0.0374363118956966,-0.0795413078365574,0.012443153156174,-0.0099366753967184,0.0379331902218671,0.0544828735173135,0.00393766167679268,-0.0531238412731904,-0.0241868905463195,-0.054697706687589,-0.0564282343426851,-0.00577347498655946,0.0454319094705544,-0.0524168368300501,-0.0151337902690907,-0.0781356664738318,0.0781029409161317,-0.0316597777379955,-0.0783601456127875,-0.0649374878336384,0.0548888318868768,-0.0307623483011324,-0.0501697162399859,-0.0421290921243509,-0.0792839760329965,0.0457219773371341,0.0699586538178654,-0.0771570660905713,-0.0238203955925165,-0.0254792333885465,0.0465186778300063,-0.0514078767464533,-0.0250205411505981,-0.0210381690883255,0.0186554479313341,-0.0670873510032368,-0.0465552956455179,0.0100968095334698,-0.0402111922764271,0.00925781380816261,0.0353991028458807,-0.0166394266191122,0.000438768731714073,0.000831012241929301,0.0109437365508376,-0.0116494515219934,0.00269534576809749,0.00904667739246351,0.0366907707400111,0.00433520015530997,0.0241091895588251,0.0115796025523821,0.0535728519286834,0.0539394733188392,0.049450510502537,0.0542888758956869,0.0196614505814675,0.0394091642458966,0.0571318098051156,0.075841054988951,-0.0660700691426499,0.0236504875606161,-0.0555668217575023,-0.011090610293248,-0.0773876816022152,0.0430886261738318,0.0263082978018108,-0.0439429772477332,-0.0268145643672042,0.0416108464270881,0.0453148366349353,-0.0714154615958293,-0.05502858026653,0.0487991753820326,0.00941555064609998,0.0359151562843077,-0.0428502761399608,-0.0678891036603083,-0.0350381663977346,0.0738404946000504,0.0164460964950016,0.0690710231610905,0.00542009715243248,-0.00998105165082079,0.0430104964054238,-0.0251293924195363,-0.0356921758296397,-0.0173280530876145,-0.0657202282481455,-0.0585603660245241,-0.0214869980986635,-0.0517902973907955,0.045090121536092,0.00294618014383417,0.0171190923159565,0.0477024400083825,-0.0339651937568398,-0.0610588690550341,-0.0762405373138564,0.0192202418759559,0.060551977372054,0.0490742993210788,0.0278047802801266,-0.0744766029689818,0.0178734747031114,-0.0427796691482792,0.0414385533153259,0.0550231985631507,-0.0306687728085875,-0.0735996131569145,0.0488636930886953,0.0657773236864141,0.0754714100041759,-0.0257162098333781,-0.0242037279644067,0.038481906335094,0.0291543977470856,0.0201040962059536,-0.0588461467525205,0.0434341695734459,0.0415437301814294,-0.000333144851184048,0.0716438721826504,0.00663385171752136,-0.0773869647818557,0.00876296449860696,-0.0256637083486019,-0.0313521585386955,0.0277040954435729,-0.0219042456624584,0.0678680833372605,0.00825607281562689,-0.0528299463413795,0.0156728635428813,0.0137794698466451,0.0450435283617319,0.0528931943946021,-0.024781976838029,0.0200667268503768,-0.0577755784884913,-0.0183815899949435,-0.0110695801354337,0.0287632150353798,-0.0229101799907676,0.043214210105694,-0.0553322911147644,0.0503980512779197,0.00254663054018589,0.0127550685092598,-0.052950058492343,0.0437145691941094,-0.0710342002618286,0.06282703905498,-0.0360649760607933,-0.0724107691610282,-0.0559578248001439,0.0446946770533429,0.0742001225027256,-0.01329639632874,0.0183118235777653,0.0185290472668265,-7.19318166710092e-05,-0.0261754899221358,0.0613510288956347,0.066272665628359,-0.0413737466006418,0.0787668615573863,-0.0385490345016816,0.06854379034999,0.00352506835177778,0.0688975491323031,0.0703360129102766,-0.00459699627226079,-0.0425797281985077,-0.0311715537082271,0.0793043083694318,-0.0179120193132721,-0.0607735025048133,0.00185093883511188,0.0748430491959877,-0.0337235609971562,-0.0344344920452845,-0.0761911511403467,-0.0508965219421762,0.00950053189392227,-0.0686019203013749,-0.02685434674232,-0.0258047910527348,-0.0744017978731551,0.0398492570034458,0.0725070325250304,0.0241272493936714,-0.040222674887731,-0.0336684574716112,0.00547827828930611,-0.0539500093338778,0.00495779592774706,0.00424513977218659,0.05758461288064,-0.0243255691716101,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.0658178208702327,-0.0563313566038065,0.0500588696310571,0.0759067577849639,0.0518925429749734,-0.0756715410927644,0.0707498069809516,-0.0618310180221829,-0.030106033063543,0.0286360995232296,-0.0209234381937065,-0.0695496883939717,-0.0708443997012658,-0.0746306596671374,0.0703086775496177,0.0636911509855143,0.0504404448580185,-0.044761634899658,-0.0213772860641485,0.0406017624403358,-0.0549978492292566,-0.0774067334073627,0.0619146629338687,0.0522798678149841,0.066663966172684,-0.0647028192806536,0.0666861274031392,0.0407490727122636,0.0358253028783134,-0.0343162461902556,0.0547388540090708,0.0588098871609242,0.0757598622868582,0.0152984704521012,0.00793581355732671,0.0608907974422401,-0.0271419077120451,-0.0165874863120669,0.0629400009582471,0.0264471505332958,-0.0679513868633431,-0.0379834373099652,0.036897462139324,-0.0587957865646089,-0.0326140969771026,0.0272061396144359,-0.0751046355790946,-0.0621736521935899,0.0624445047147779,-0.0164819216432431,0.0584281103212517,0.0647932145821144,-0.0138886550506058,0.0403427731806146,0.0597265232259997,-0.0272246888779219,0.0556399539744667,0.0464126505546331,-0.0664756162401641,0.0114652567782743,-0.0679035956356226,0.0682632378434126,-0.00972485613530728,-0.0721437334232702,0.00356170822100793,0.0782109574220194,0.0687470640189699,-0.0217880935695898,0.0640296162590522,0.0516870649027112,0.00286695096775281,-0.00632791582789641,-0.0662963724817598,-0.040235586967429,0.0148762976820005,-0.0189104693843566,0.066970552647007,0.0321053442408369,-0.00108412157794652,0.0494150573617849,-0.0767102596148431,-0.0226560113312006,0.0568617127728005,-0.0105989146654489,-0.062313238150586,0.0365882359988001,0.0421763964566292,0.0733267158238807,0.00300088647892739,0.055700780290971,0.00479197252764925,0.0150972908433048,0.0439640180598777,-0.054516315969879,0.022618413391811,-0.0572728107018735,-0.0244032902011663,-0.0383109099223795,0.0168266707364594,0.0173158894187379,-0.0224132672429147,0.0674314556212311,0.0146710545265447,0.00483120601849221,-0.0061529317899388,0.0509080148352813,-0.0575229979201793,0.0563938426861511,0.0680953267906305,0.077270823864858,-0.0658073097401333,-0.0690422282316919,-0.0472917724248449,0.0226031510264627,0.0329475054670812,-0.0338470014714855,0.00347060505462373,0.0484186518045229,-0.0534132137025768,-9.4295702918572e-05,0.00535699721675226,-0.035495290139455,-0.0313867033791666,0.0633975579512294,-0.010011606109334,0.0712317100871502,-0.0738752528251499,0.0455851036894997,-0.0470791999097351,0.0229514179113095;
//    springCoefb<<-3.86450119769613,3.67396148701663,-4.96849191583509,4.33734622681848,4.17603743289614,-1.33568579837623,1.15012467871814,3.95741490354904,-5.68490169464001,1.60881794182029,-0.876656403518422,0.259682627829472,-3.40747995045241,-2.18698379600172,-3.72341609895633,-3.57067374939865,-5.12424297288473,-3.18055398600025,0.417737634523928,4.8235314827992,-2.5166620583604,2.33178977414935,-4.53172942699384,0.906239655311621,4.61007161467486,6.16366426825652,-2.87411114577893,4.65809081727067,-3.68895724557279,0.0398834617873903,0.620399915144564,5.1743609479423,-2.25685082540471,1.93509319930948,3.22852197476078,-4.36399859250857,-5.68377780491843,-1.90453855237275,5.87660151104047,-5.08549429955844,5.98746458944755,-1.28324009247795,1.45737352827104,-3.70320056685654,2.81296131152033,4.01714263516638,-0.990689116255186,3.9719035386356,-1.76097794466199,-0.931188792323655,4.7886848273943,4.8221907184018,-5.02425520974346,-2.26725924159576,0.546693422196206,-0.0183688877131311,-1.37746925424591,-3.18839334253656,4.88263154649922,4.586231539247,4.27729203625966,0.626629874566872,2.27658908341489,-4.64607910140807,2.86791570537105,-5.18858848726607,4.55014480331692,2.9092112956244,-2.31733057297406,-2.29792433884564,-3.71580874952019,3.29165047590883,-5.62185617605571,1.51884062849406,1.07067794109222,-5.74215950563337,0.708732434591533,5.82309725519261,-2.13598771901082,5.23093968992955,-1.39127674298272,-3.63048809746819,3.76994520833135,0.832505528858159,0.385437866787722,-1.96654657532412,6.13246956541236,5.29115381254182,1.12824528213933,4.73191590605991,3.59420015178882,-0.877647887452685,-0.924639425224893,-0.412395970647958,0.759458216990917,-4.33990892570552,0.682200742085968,-0.973582185543834,4.85248756991888,4.6480553691119,3.01167867561052,-5.14650637960131,1.65652063916906,3.67300769955481,2.65551944889275,-3.55598661973871,4.21403339977311,-2.91893332236739,-4.01607457039777,-4.20513952508938,-3.97117883243784,0.87583388661951,-1.55244242255757,6.08195157589351,-5.53969807195834,5.11618064423015,-2.16778020528228,4.55860292152544,0.706917427270735,-0.527725255953519,-0.322841148056009,5.29013766380914,5.94964348008737,2.32649434433559,6.00465161784824,-5.52977231924701,-4.56082833002572,2.85335845942771,0.988860129289455,-2.24536976258571,3.84903450826306,-2.74962286407637,1.82795550024617,1.6256001749435,-3.44664833496869,-0.787624398580966,5.25927800430962,-3.37926496591605,3.72316823967611,2.58035980433033,2.54019413687263,3.00822696142147,5.43975856755787,2.60790289401073,-0.295271064611664,-4.88084775784971,-4.20844502939635,4.02247423416613,5.96094036367573,2.78165760372606,-2.78843622763906,-0.645085984380277,1.78861006168353,6.14382046948618,-4.60177684589636,1.51007647953178,-5.66913705561251,-2.87941997592207,5.98029692774147,1.60290827367695,1.15839546734389,-4.35390095862913,5.13647060960058,-3.29683423826162,3.55488441631438,-4.59336292536811,2.19872656315742,2.530977220624,-1.68944269128416,-0.361290403018153,-1.17184818089734,-5.4324337602632,-3.63624864744836,-2.01527481919115,3.45865433374752,2.35166548793998,-0.612937377040858,5.53339451020285,0.0909545162544419,-0.935182213365127,2.03186690807724,3.58570349446706,5.95180612610385,-2.4627082360909,3.44633875810156,-6.18204504949343,5.33055344344088,4.06038690834072,-2.77827981956384,-2.87238202345131,-0.619890023834005,4.66330084778005,-0.943097782080431,-1.76660461423343,-4.91671859633324,-3.67139857161773,-0.076782339601541,3.56519316682417,5.14276384900627,4.5169601691143,-3.07928244204565,-2.31226953774274,5.36771161470277,-0.432345889494007,1.95564084306611,2.54318074259862,-4.3638656074057,-4.94048173397475,1.7933900469498,2.01027410884874,0.407521252660124,-2.45792825082464,-0.687207602535872,-1.17274695093695,1.36254871893614,-3.52405404443431,-1.07160680043039,0.409916956525344,-5.74685234194526,2.43329858000578,3.82072013307404,-0.0835571657792664,0.813414221934155,-3.40556285485807,4.4330234199873,2.17988082560091,-0.79377622062412,-1.92694412546591,-0.538111213426586,-1.93419757746952,-3.69316915635161,2.66579154452776,2.0367180906394,-4.60864274750051,-4.04973955081792,-2.29082627214616,4.21772319509811,-2.13041995822362,-0.948122806120907,-0.272071957952097,6.1630393564768,5.74258364653922,3.55318499122327,-0.807353451910744,-1.71334851024941,-1.36745149569227,1.95177770365494,3.49822989517188,5.32565066668475,2.48811056170968,-0.351656730674014,2.86318559390712,-3.87863180406958,-5.82142770873986,5.74080794490072,-5.72879358993395,2.64163831686106,-1.33615348157507,5.62309623073544,-4.17965809656553,3.01283414095541,-3.09774938915601,4.09833157945901,-4.45983733279237,4.38703961039144,6.222876169116,5.39500339261633,-2.52784178346676,2.58829038843656,4.29048320347157,-3.55891374906826,6.15197924540229,-1.36853877735244,5.51428064317356,0.833913712280088,-2.52541953056468,2.13240475181388,6.06555108245924,-2.58978894902704,-0.581028044578514,3.28479813262921,3.59677566583801,-0.755150888844877,-3.03782608294973,0.795958958730565,0.997268549245814,-2.59027350259689,-4.79681275329447,4.1383073133398,1.06370354363343,-1.01475946668824,5.37583775584291,4.24913935447742,-3.19961308722923,-5.36718477694946,2.35299375901719,-3.25992212396491,-6.2553665901848,6.1083371814021,5.61155347032332,4.31830181328677,-3.73376177351783,5.48034750987394,-3.33342216991734,-4.50266633034427,0.0310760163023503,0.424343499517975,3.9129236214696,-0.186558101238407,4.11773975049093,-2.95128962896058,-3.18494517446087,1.43133021047728,2.57674468804621,5.20067858627806,-6.1831742934951,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,-1.86624674250335,3.61694724017608,3.16654841230712,-3.50020186774831,5.78247426784757,-4.84539688217807,5.13597709712056,-3.76063306196901,-4.81757827236287,-5.3107692849138,2.89459043030069,-3.82998730933537,-1.27757149567505,-1.93629822246639,2.86057870319148,4.33289757493666,3.7212432150659,-2.7154157788901,2.02822266216881,5.25053533098438,-2.05733900903557,3.71929467258504,-1.29901604863953,4.00948049067158,3.49566395445844,-3.87413420294964,-0.814363134953018,-1.13390855392051,-0.132305036154921,4.07865340641918,-1.20756497996135,-5.67315378971611,3.60556657449736,-2.16231139627231,1.6799593601027,-6.03598319715729,-1.45293454920339,-1.23471358501071,-6.09307128039791,-4.58914961603963,-0.26229766407284,3.08470434990278,-2.13595171952332,4.74331604025211,-5.13477907841529,-5.55855822218352,2.79302841518877,4.86964934250228,-1.99078880107362,-1.4619341284941,3.83699946763499,2.90383938221128,-4.02582465590906,-3.74520178685621,-0.0386473194376096,5.75302449854938,-1.33615078980585,5.43017474560937,-1.6640692612228,4.81472937403923,3.22564294617688,3.41155095881585,5.4247607901748,0.548024314822573,-5.03394563745646,0.821534944425819,0.795226317665279,0.0300845206647421,3.96434777754175,0.985340237267368,1.49034060315219,-0.6754763105094,-2.21314061868153,5.63757408362887,-2.21534547025729,-1.06473449124515,-6.20416934440633,4.70791372758925,-2.47827035459454,-1.91177294547994,-0.884251189414299,-4.92445609281122,-5.9599007614408,1.37310935467664,-2.38647267966743,0.284637119121588,0.842948653226024,2.56056173052672,-0.568373341120713,5.46206459200322,1.09210590456595,-3.6259156007955,2.59043035081907,1.3940974479496,1.14890059666481,5.02898847525872,-5.93975277519345,0.208744968082346,1.64427444680061,-6.04664379701433,-2.72603572718022,3.47833679606309,2.05592782561468,-0.736004088604331,-5.48585112463027,2.08016178295133,-6.07073525442908,5.49252503872056,-0.897913872670351,4.8649391267239,-0.386430778530009,-3.67537460148936,2.32622207187887,1.54348373475234,6.20180004328616,-0.911167348488062,1.04370269768405,-5.34840537980092,5.23027941649461,-2.35946564997576,1.30411456968309,-0.795897492765567,3.33436945617133,6.01092648885484,-5.68498525066764,-1.79991515301554,4.75672975826189,-5.34155282586109,4.69201501506681,0.117819005062502;


//    generation1112
//    best = 15.4793
//    springCoefa<<-0.0787365553522187,0.0397860520993408,0.0285811549930746,0.0148356401616873,0.0225603931014757,0.00942131270161891,0.0190325777693803,0.0233728654418946,0.0583954952556619,-0.0483040296883807,0.0242124729536801,0.00489416932481162,0.0665881243375549,0.0251908738910303,0.0741066694232201,0.0670986463979148,-0.0345643036600967,-0.0259200940029323,-0.0768512000128865,-0.066696290442113,0.0769517727367604,0.0440311868354987,0.00310712193810026,0.078899163230741,0.0788787391934348,0.0771249259249889,0.0117853437113109,0.0767958074048142,0.0437036912719271,0.0292698902400536,0.0590928833160164,0.0449671359942142,-0.0638690479769693,0.0677258177135295,0.0741325108270997,0.0139300205356702,-0.0556104941552554,0.078835353850776,0.00860470357830995,-0.0199256604630154,0.0330242166540873,0.0383049539282475,0.0671731600105638,0.0100578547688471,-0.0763951291127108,0.0609859394193561,0.0691503421495137,0.0247812721674833,0.0637547364010265,-0.0723889178002202,-0.0519388390946849,0.0624782460846371,-0.0431796303219998,0.00900752563448974,-0.0400732179545207,0.0376800872188434,0.021020399800878,0.0117542496590411,0.0336049248422409,0.0625928757258658,0.0278127577633072,0.0702474944993143,0.0435802303085012,-0.0415860856145556,0.0662081225558886,-0.056616993684609,0.0702688369854674,0.0669085326922944,-0.0577816398338329,0.00498168028038696,0.064275964950129,0.0082106890540437,-0.0675601910741814,-0.0297229902212149,0.0547141934895084,-0.0238513961172902,0.0504024790089589,-0.0269918353608771,0.0652930952540101,0.0246078541616946,0.0783054675153948,0.0217986215452878,-0.0539794790623614,-0.0749792838073239,0.0455257733005685,0.0198429829005026,0.0347960266260412,-0.0795648236570716,-0.060654577007822,-0.0186328556102854,0.079048270806227,0.0673910252187208,-0.0283853611854768,0.0768474262419641,4.22997400361602e-05,-0.0458803351064587,0.0660115074301192,-0.00968886334900225,0.047408684808348,-0.0717701324782195,0.035807761529371,-0.0174382546997807,0.0189957893067148,0.0482475704551896,0.0328387551535101,0.0190584340176443,-0.0556038256621006,0.00324123408796323,0.0127268029229173,0.0472632276808765,-0.052150911824848,0.020264137449631,0.0174305656260953,-0.0261303908872094,0.035257932212356,-0.017043661147842,0.0726981328223454,-0.0570304115009636,0.0253911636689594,0.0197565364348446,0.016967082249871,0.0102833213978636,-0.0188242315961161,0.0559513717032743,-0.0549317130329701,0.0612180681439201,0.0193526401982814,0.0226683045564289,0.0065515928256988,0.00290657720210824,-0.0606903381555762,-0.0726630338247228,-0.0690636828490829,0.0193897279593573,0.0555845366304668,0.0437750723044272,0.0468450873002075,0.0394696528074704,-0.0329836936821154,0.0769758056090101,-0.0703300195887359,-0.00513460550696338,0.0571396772193301,0.0142195497380376,0.0164515214097471,0.0308472113656125,0.0694353171229428,0.0780693864068339,0.0725435366842717,0.00753803458415812,-0.0130345523790617,0.0691651284023265,-0.0578139640287561,-0.0610438296855631,-0.0769839821741842,-0.00349377188994259,-0.0491218169634798,0.0680368567239194,0.0104728014513025,0.0277090399468825,0.0643284656034449,0.0148744277499296,0.0319480187967178,0.0469556046268861,0.0225398461253102,0.0540060857562377,-0.0767194037869197,-0.0488786590699473,0.0358852142383299,0.0235482812711195,0.00388935815724048,-0.0344991294827774,0.0072645609730069,0.0182211471806379,-0.0330087505062151,0.066253237997802,-0.007020162309995,0.0507489917067276,0.0579286789791326,-0.0554916582887488,-0.0443672990633023,-0.0351058734744349,0.0638490229759595,0.062662030843652,-0.016149703159998,0.0502813330899371,0.054324965017999,0.0552459516260987,0.0749330744627788,0.075759350480195,0.058129147798243,0.0282569440958355,-0.0268196486434059,-0.022496843981788,0.0128301607010894,0.00652780320799341,0.00679419323244019,0.0513311464426714,-0.0737902242661408,0.00482908282653852,-0.028651634505322,-0.00741968561309375,-0.0113375247695192,-0.0309082336495203,-0.0691985385069617,0.0356537247242656,-0.0510489410772216,0.00378129925754912,-0.0162516089231947,-0.0731202621725948,0.0282896409688004,0.019381092013503,0.00626741023776652,0.0410187031411662,0.0771998289214446,0.0356241612674781,0.00990672205602387,-0.0794443978366649,0.0199252965021531,-0.0576569127932456,0.00239666443429732,0.000663004027198092,-0.0266708465417246,0.0555770157908914,0.0329926150670633,0.032585694171761,-0.017895181075621,0.070577879464054,0.0374795822822829,-0.0116854053417618,-0.00815244164697474,0.00814443644515446,0.0608949091196502,0.0605100336580118,0.05723620287014,0.0716963706126885,0.0161637583077716,-0.0738127382815875,-0.00452233020426814,0.0799121494590827,-0.0669330003796765,0.0132565554729351,0.0192932413980799,-0.0151591360267061,0.0587106682500125,0.0164930703195245,-0.0595349748337339,0.0651585560967953,-0.0119821358900434,0.0313351379853371,0.0102307653847294,0.0735335526123895,-0.0783080907903184,0.00670429681703182,0.0459915443351453,0.029551836526744,0.01614561294026,-0.0519036367404757,0.0165703120532307,0.0635226350869781,0.0793850834321845,-0.0715821296682498,0.00138656943146689,0.0789834328713299,0.068927903989762,-0.0416776768684749,0.0791554930915526,0.00509166229753367,-0.0354904150755566,0.0487799919627843,0.0385185939230934,0.0429319663492875,-0.0717527819386464,-0.0557029469943153,0.0424174485180608,0.0275694849097959,0.0283013731032796,0.0628824737588328,0.065324034305884,0.00756486963832979,-0.0580133232371944,0.000926739197192052,0.040198942275384,0.0671719744676738,-0.0315740670038266,0.00459441494531404,0.0406146016064966,-0.0549752013268765,0.0632378235381273,0.0168862972014782,0.0632471188079785,0.0675372461255414,0.00629654168699845,-0.0716245620658736,0.0733686844228621,-0.0731565964935145,0.0207637777727245,-0.063657464228411,0.0332426104476874,0.0721893606780537,-0.059819667255422,0.0611593644201538,0.0371006349274426,0.008756611388529,-0.0662300293828501,0.0403580811924986,0.0173801001790549,-0.0183408645998411,0.0773245660948216,0.0381073966427275,-0.0354915834756063,0.0505520611814888,0.0779140615453543,0.0467182230794422,0.00691535268300928,-0.00752095301054463,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.072972533086243,-0.0563313566038065,0.00493563063820985,0.0168609884832338,-0.0755998318621888,-0.0764228115959199,0.0651290465822113,0.0677687595727708,-0.0255279210701249,0.0154125595537073,0.0310159908752032,0.0324546150026772,0.0672664034017367,-0.014407729028914,0.0427632872586899,-0.0203874216696189,-0.0329939846475581,0.0582008473473605,0.049530052081463,0.0420176874102258,-0.0526738212130376,0.0566280047414426,-0.0299665013467737,0.0202939150204388,-0.00878208395502627,0.0613053419223006,0.0651859232358668,0.0415810448984855,0.0684167652523223,0.0225791594861915,-0.0764551426453773,0.0119980309572516,0.0657973458551789,-0.0340542982304722,-0.0545746671290019,-0.0363075187598856,-0.0704504885480043,0.0251474474115052,0.0239696648456015,0.00888312433328625,-0.0184900807302865,-0.026193907384851,0.0346668059912821,-0.0502069033543611,0.0727939190186036,-0.0652426554193919,-0.073405683857112,0.058832445019454,0.0433435391706107,-0.0276349627727805,-0.0627231468924895,0.0216381766375332,0.0484023522030692,0.0511634874649726,0.0491580558268187,-0.0303264040268615,0.0372968669781912,0.0201023868931081,0.0588998787006828,0.0257136321560077,-0.0356257363947228,0.0141284905145072,-0.0589804615727535,-0.0498283906140497,-0.0516095621751666,-0.0335551287017554,-0.00613590937393527,-0.0420600507231709,0.0715923187097499,-0.0621662446028396,-0.0562829418369955,-0.0100930846902044,-0.00836015198769055,0.00148987960139749,0.00289533462510226,0.0398096824344578,0.0162472241820056,0.00948965076799023,-0.0518638096618763,0.0600174140898688,0.0472304060902588,-0.0355273761393164,0.0202883375903072,0.0234392278378081,0.0422216261339447,-0.0518900550770993,0.0370514884947964,0.0149342752690121,0.0455472842443489,0.0482517623753528,0.0192839286938701,0.05911724067252,-0.00100677106576356,-0.0409349269237998,0.0659956001425141,0.0674249651969527,-0.0747426529204206,0.0431735729416369,0.0108531262403601,0.00681786573268896,4.67374548533606e-05,0.0628869602749529,-0.0176661403186927,0.0786480859660767,0.0524744645035985,-0.0418882346720846,0.00538064945739723,0.0144900607737209,0.0412227917798612,0.0745418645046977,0.0701532715513153,0.000597512655533077,-0.0725604402984308,0.0226924941421917,0.0577535415337205,0.0666303270992964,0.0287497000157599,-0.0541356280325612,0.0614123152109852,-0.0198186871129175,0.0215594750847957,0.0270599850426697,0.0423853336672181,-0.0427899707680521,0.0499074528313328,0.0793047787834378,0.0473748562227119,0.0524015135663986,-0.0379093903107147,0.0793937288558049;
//    springCoefb<<-2.43192959884493,6.08514030121006,4.04035678062355,4.0939778826623,0.885759964104128,-4.40632201143234,-3.84532036998163,-5.33669165449445,-0.201902543387136,-5.74421405700816,1.84616840882136,3.98087731087271,1.11369296464456,0.934977743915378,2.02614908115595,0.959073349058476,5.15053093739066,3.006446367785,-1.6732576622727,-1.23768408693475,1.27684178065925,6.04271727318241,3.77882575946868,-6.08085982355576,0.652931759876538,-5.48543869292053,1.81651248891278,5.96869830004477,1.31094790673134,3.23199520554859,0.163817214486977,5.16220350788641,-2.7309598098749,5.48964219698341,5.39277565112383,3.98484668574114,1.44712120543646,5.41086102056708,5.39994485533409,0.946468807561242,4.0832679471894,0.279088417230802,-4.91215907281787,4.6082018532571,0.92838679515286,-1.64664393517429,3.42329263641773,5.70322596816668,-2.65748894621183,1.29342632652648,0.568687125340165,-2.43705856028738,-1.92000280751018,-1.8801132043928,-0.517201690556528,1.68458566290145,4.00263656362903,2.50000791231667,0.452871199207279,-4.12575887624951,1.3802421085177,2.0898833298572,1.09336188408456,4.51445899470899,0.99940722926591,-2.21682713707008,1.5502626395456,5.9937160473333,-3.08915132235467,5.79470203052607,5.86985869733578,1.12879047483547,2.60477802497352,0.454742152576666,0.517470415219057,-5.5002405095043,-0.520598258586578,-5.43747913056799,-2.17154290659146,-4.77267751366817,0.0886738094516904,2.80624253330262,6.11985611039204,-5.33693277264827,-6.2580654340777,1.58892567487088,1.83628608125866,-5.35466976461089,2.3785831099318,-3.87135237064131,3.14285901158726,4.07157566778108,4.50171615921588,4.25724858523384,-4.71109470247867,-0.420077554397874,2.01939375860174,3.12235313706693,3.04526415542555,5.21342763624706,-3.62827223889698,-0.737268266315717,1.96095282746985,5.25969098607654,6.00065908626095,3.66813714062273,6.04263567657224,-0.803124378177302,1.71191742881433,0.288838383269688,0.707383314006201,3.10460104757031,-2.49262183233593,0.54405421854657,4.53851627650281,-2.46750206641364,2.16747008949913,-0.845862873766654,2.96766701262683,1.696825034699,1.26885943120368,1.67349275491913,-1.65457166366446,-0.215499085192086,-0.443091929617129,-0.0824811602914577,3.72336076311685,5.3500157115321,4.86531505313461,1.28618349288566,-5.77664094062,-0.588400267973177,0.716840377031068,4.52819929759633,-1.61189448189663,0.434314257440346,4.15175172698769,2.17904391568207,5.91437508511471,0.300796283741831,2.01864875880472,0.338573193269244,5.38597896424107,5.14448955286199,3.82100555758673,4.50086791031566,3.6169715871486,3.34167081581541,5.25594889042253,3.6242449893384,2.869005561767,5.92319047857269,4.68127868871957,-3.14173898982101,-2.97035686401395,-0.230860614223697,-0.657414860877198,0.394065895226855,1.19275490160846,5.56181821803563,-2.91342920436725,1.82579686459939,-1.63544369410512,1.79662060059091,3.02465413262939,-1.00109808835509,4.72917326240044,4.64586125535466,2.46719717703016,5.95293752900116,-0.728069876461334,5.10586141337541,2.35273521470591,1.07370842850819,-1.93135481169837,4.37060868490869,-5.19487679371209,0.391702183144271,-1.21653276126522,-6.13711309553413,-3.61258757467371,-4.63071240534989,1.85101678249404,0.308870603002905,-1.4892661951709,5.49822140170814,0.837830493970494,5.16937915605619,5.61545234316673,3.98285231687726,4.89567875026041,-6.23992361757923,-5.34797394185981,-2.9636930323549,1.39715120233631,-2.12644323785466,3.40947531887144,1.80457371889052,4.09084134539069,-4.52002808382579,1.6512736377946,5.63214168372345,3.85473970142402,4.58615899606649,0.42266490637997,-4.35980031027435,0.0279554131341941,1.51097331266788,0.96973753156521,5.09460785186898,1.65704541713375,3.6403351568915,1.05442214408828,0.661927881806177,4.70902626508569,-1.02537075450348,0.192632458299389,-2.51630261939429,1.84969643455691,3.83759113021315,1.5746427306926,5.27170062187801,-4.70822915029589,2.50985398883279,2.6840672558054,-6.16522383252989,-5.89977445487354,3.05056480308843,2.21451166528256,4.47425209051715,5.5492140964325,3.75304063429704,3.82320856838893,3.12076859200485,2.05601442451186,-2.03731173108278,5.0441534817305,-4.19921536235395,5.75684678158511,-0.269294192555964,-5.38779271633664,5.43889699630295,-2.91214424151613,-4.92389726981756,3.36810651670152,-4.48630318228211,0.33391717567896,-4.53575025726799,1.0605535709224,4.47792557723166,-3.00405579583632,3.00532537434112,-0.591793793133674,2.71423009991725,5.14505028459612,-2.34202100170558,-1.31113837281043,5.52846103557426,-4.83914922729897,3.79489135958476,3.12888332563001,-5.57312033671814,1.83308250263569,3.19973161421444,3.83083346113838,-5.96668672335076,3.33493105350654,2.5918017370172,-3.88271688570471,3.38475512344542,4.01203987647953,2.89624094641009,-0.173391873853497,-0.589636697054897,-1.62803646600724,-1.48343160398213,0.848696330955331,4.98906590967172,4.47871139269655,1.19791773758934,-4.04150873496511,3.90776009537591,3.43991137256052,2.61626586932033,1.49626440893562,2.84893529190212,4.7158110969248,-3.53724639953568,4.50491096743946,1.49275078504112,-5.77939771606788,4.17997546079867,4.1568390815018,-6.27213179494185,-5.02828562505671,-3.57958680001224,3.34456559070511,3.46028085870612,0.799359763371303,1.47706369856494,2.56591515626386,3.63946459053563,5.86683490040684,5.73396412750752,4.27641316775727,4.53509332871182,5.61569700322845,2.9393886033297,-2.31855975768557,-5.79668449207325,-5.13523515189644,4.72865749192539,-1.49676813258366,4.1860763159865,3.37062530542857,1.21079109463131,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,5.12690620844648,3.61694724017608,3.05138930367418,-4.01275011252449,-2.72372981137355,-0.0180124155448347,4.12303341309968,-3.74873780023526,-0.887138111904024,-5.41990213502721,-5.9774984601979,0.552717896425435,0.921612258771997,3.17838284590814,1.86797359438216,-1.47430087116196,1.35138812294626,-3.90908982234295,0.382622065522161,3.58663702089682,-1.89477216410266,0.842265657489956,5.32264774649661,4.07279963344888,5.93452553491607,4.58997732369303,0.786062578409914,0.325200986907274,-5.9917002755706,0.588543628002016,-3.55057119553441,5.24768262746627,-1.25341481552174,-2.18636580089857,1.02549154472145,-2.97758908615075,-3.02601767758674,-0.79860588073586,-1.31513570862061,-5.10264541624672,-4.04774808632133,4.46998918241493,1.57982980509186,5.8374249517381,0.290319230405863,1.86919610233567,0.185816115052586,2.65571266392906,3.43118334743,1.54324804109866,-0.19298985557722,-5.47214491078457,2.00854924420448,2.55367598272938,2.67066286720342,-4.40899082454341,-4.04429312445419,5.6046481394415,-1.30713353515645,-3.75280819417312,5.43198984070525,0.316535709131635,3.81446074559214,-2.10461017481648,5.52229986841057,-1.44323290968642,1.20098593903276,-3.78690300917618,4.04134640957772,6.16903543626382,-0.65473740484626,-5.61395256815843,4.35583941282708,5.25665177966897,5.8310232749945,1.77763132207111,0.842662682004633,-0.266345809952913,1.76766768752563,2.41652541542875,4.29474875067796,5.02212468610824,-3.23113956818491,-4.46801644912111,3.88677049527466,-5.81971337358793,2.74959988627578,-3.37151591823711,2.76592180759093,-3.16280122963266,2.80490867323494,2.55004818904424,0.689441449147486,-0.0351284046447834,6.06388826542954,3.55993961157682,0.535486416593954,3.05170925752847,-1.00278238102184,5.46116454437677,5.81344950372511,-1.74189744386769,-0.332778030647898,-4.03931177144349,4.25502161714161,4.25874652459808,1.4159745097112,3.55494616590075,4.0525786915951,5.09901468640174,4.31990406550899,1.04140371424045,5.82363717143333,-2.31027463831456,-0.906303222026946,-0.757257994087879,-6.03794887918608,0.696102823104212,1.00660345861309,4.57369870271486,5.17365517040888,0.38104645055342,3.68186372776027,3.74661889802106,3.58529704079631,4.04925640652263,4.59815072532965,1.95431786572231,-5.33829054200515,1.63610028456171;

//
//    generation1199
//    best = 15.5074
    springCoefa<<-0.0787365553522187,0.0397860520993408,0.0285811549930746,0.0148356401616873,0.0225603931014757,0.00942131270161891,0.0190325777693803,0.0233728654418946,0.0583954952556619,-0.0483040296883807,0.0242124729536801,0.00489416932481162,0.0665881243375549,0.0251908738910303,0.0741066694232201,0.0670986463979148,-0.0345643036600967,-0.0259200940029323,-0.0768512000128865,-0.066696290442113,0.0769517727367604,0.0440311868354987,0.00310712193810026,0.078899163230741,0.0788787391934348,0.0771249259249889,0.0117853437113109,0.0767958074048142,0.0437036912719271,0.0292698902400536,0.0590928833160164,0.0449671359942142,-0.0638690479769693,0.0677258177135295,0.0741325108270997,0.0139300205356702,-0.0556104941552554,0.078835353850776,0.00860470357830995,-0.0199256604630154,0.0330242166540873,0.0383049539282475,0.0671731600105638,0.0100578547688471,-0.0763951291127108,0.0609859394193561,0.0691503421495137,0.0247812721674833,0.0637547364010265,-0.0723889178002202,-0.0519388390946849,0.0624782460846371,-0.0431796303219998,0.00900752563448974,-0.0400732179545207,0.0376800872188434,0.021020399800878,0.0117542496590411,0.0336049248422409,0.0625928757258658,0.0278127577633072,0.0702474944993143,0.0435802303085012,-0.0415860856145556,0.0662081225558886,-0.056616993684609,0.0702688369854674,0.0669085326922944,-0.0577816398338329,0.00498168028038696,0.064275964950129,0.0082106890540437,-0.0675601910741814,-0.0297229902212149,0.0547141934895084,-0.0238513961172902,0.0504024790089589,-0.0269918353608771,0.0634395048810274,0.0246078541616946,0.0783054675153948,0.0217986215452878,-0.0539794790623614,-0.0749792838073239,0.0455257733005685,0.0198429829005026,0.0347960266260412,-0.0795648236570716,-0.060654577007822,-0.0186328556102854,0.079048270806227,0.0673910252187208,-0.0283853611854768,0.0768474262419641,4.22997400361602e-05,-0.0458803351064587,0.0660115074301192,-0.00968886334900225,0.047408684808348,-0.0717701324782195,0.035807761529371,-0.0174382546997807,0.0189957893067148,0.0482475704551896,0.0328387551535101,0.0190584340176443,-0.0556038256621006,0.00324123408796323,0.0127268029229173,0.0472632276808765,-0.052150911824848,0.0281599258042537,0.0174305656260953,0.0441711298954484,0.035257932212356,-0.017043661147842,0.0726981328223454,-0.0570304115009636,0.0253911636689594,0.0197565364348446,0.016967082249871,0.0102833213978636,-0.0188242315961161,0.0559513717032743,-0.0549317130329701,0.0612180681439201,0.0193526401982814,-0.0689202056773567,0.0065515928256988,0.00290657720210824,-0.0606903381555762,-0.0726630338247228,-0.0690636828490829,0.0193897279593573,0.0555845366304668,0.0437750723044272,0.0468450873002075,0.0394696528074704,-0.0329836936821154,0.0769758056090101,-0.0703300195887359,-0.00513460550696338,0.0571396772193301,0.0142195497380376,0.0164515214097471,0.0308472113656125,0.0694353171229428,0.0780693864068339,0.0725435366842717,0.00753803458415812,0.0422359183199151,0.0691651284023265,-0.0578139640287561,-0.0610438296855631,-0.0769839821741842,-0.00349377188994259,-0.0491218169634798,0.0680368567239194,0.0104728014513025,0.0277090399468825,0.0643284656034449,0.0148744277499296,0.0319480187967178,0.0469556046268861,0.0225398461253102,0.0540060857562377,-0.0767194037869197,-0.0488786590699473,0.0358852142383299,0.0235482812711195,0.00388935815724048,-0.0344991294827774,0.0072645609730069,0.0182211471806379,-0.0330087505062151,0.066253237997802,-0.007020162309995,0.0507489917067276,0.0579286789791326,-0.0554916582887488,-0.0443672990633023,-0.0351058734744349,0.0638490229759595,0.062662030843652,-0.016149703159998,0.0502813330899371,0.054324965017999,0.0552459516260987,0.0749330744627788,0.075759350480195,0.058129147798243,0.0282569440958355,-0.0268196486434059,-0.022496843981788,0.0128301607010894,0.00652780320799341,0.000521415521053668,0.0513311464426714,-0.0737902242661408,0.00482908282653852,-0.028651634505322,-0.00741968561309375,-0.0113375247695192,-0.0309082336495203,-0.0691985385069617,0.0356537247242656,-0.0510489410772216,0.00378129925754912,-0.0162516089231947,-0.0731202621725948,0.0282896409688004,0.019381092013503,0.00626741023776652,0.0410187031411662,0.0771998289214446,0.0356241612674781,0.00990672205602387,-0.0794443978366649,0.0199252965021531,-0.0576569127932456,0.00239666443429732,0.000663004027198092,-0.0266708465417246,0.0555770157908914,0.0329926150670633,0.032585694171761,-0.017895181075621,0.070577879464054,0.0374795822822829,-0.0116854053417618,-0.00815244164697474,0.00814443644515446,0.0608949091196502,0.0605100336580118,0.05723620287014,0.0716963706126885,0.0161637583077716,-0.0738127382815875,-0.00452233020426814,0.0799121494590827,-0.0669330003796765,0.00364985590321194,0.0192932413980799,-0.0151591360267061,0.0587106682500125,0.0164930703195245,-0.0595349748337339,0.0651585560967953,-0.0119821358900434,0.0313351379853371,0.0102307653847294,0.0735335526123895,-0.0783080907903184,0.00670429681703182,0.0459915443351453,0.029551836526744,0.01614561294026,-0.0519036367404757,0.0165703120532307,0.0635226350869781,0.0793850834321845,-0.0715821296682498,0.00138656943146689,0.0789834328713299,0.068927903989762,-0.0416776768684749,0.0791554930915526,0.00509166229753367,-0.0354904150755566,0.0487799919627843,0.0385185939230934,0.0429319663492875,-0.0717527819386464,-0.0557029469943153,0.0424174485180608,0.0275694849097959,0.0001801841303802,0.0628824737588328,0.065324034305884,0.00756486963832979,-0.0580133232371944,0.000926739197192052,0.040198942275384,0.0671719744676738,-0.0315740670038266,0.00459441494531404,0.0406146016064966,-0.0549752013268765,0.0632378235381273,0.0168862972014782,0.0632471188079785,0.0675372461255414,0.00629654168699845,-0.0716245620658736,0.0733686844228621,-0.0731565964935145,0.0207637777727245,0.0562759769744742,0.0332426104476874,0.0721893606780537,-0.059819667255422,0.0611593644201538,0.0371006349274426,0.008756611388529,-0.0662300293828501,0.0403580811924986,0.0173801001790549,-0.0183408645998411,0.0773245660948216,0.0381073966427275,-0.0354915834756063,0.0505520611814888,0.0779140615453543,0.0467182230794422,0.00691535268300928,-0.00752095301054463,-0.0322656234876558,-0.0783907326303379,0.0315614002717479,-0.0782915799311789,0.072972533086243,-0.0563313566038065,0.00493563063820985,0.0168609884832338,-0.0755998318621888,-0.0764228115959199,0.0651290465822113,0.0677687595727708,-0.0255279210701249,0.0154125595537073,0.0310159908752032,0.0324546150026772,0.0672664034017367,-0.014407729028914,0.0427632872586899,-0.0203874216696189,-0.0329939846475581,0.0582008473473605,0.049530052081463,0.0420176874102258,-0.0526738212130376,0.0566280047414426,-0.0299665013467737,0.0202939150204388,-0.00878208395502627,0.0613053419223006,0.0651859232358668,0.0415810448984855,0.0684167652523223,0.0225791594861915,-0.0764551426453773,-0.00469409365425543,0.0657973458551789,-0.0340542982304722,-0.0545746671290019,-0.0363075187598856,-0.0704504885480043,0.0251474474115052,0.0239696648456015,0.00888312433328625,-0.0184900807302865,-0.026193907384851,0.0346668059912821,-0.0502069033543611,0.0727939190186036,-0.0652426554193919,-0.073405683857112,0.058832445019454,0.0433435391706107,-0.0276349627727805,-0.0627231468924895,0.0216381766375332,0.0484023522030692,0.0511634874649726,0.0491580558268187,-0.0303264040268615,0.0372968669781912,0.0201023868931081,0.0588998787006828,0.0257136321560077,-0.0356257363947228,0.0141284905145072,-0.0589804615727535,-0.0498283906140497,-0.0516095621751666,-0.0335551287017554,-0.00613590937393527,-0.0420600507231709,0.0715923187097499,-0.0621662446028396,-0.0562829418369955,-0.0100930846902044,0.0212442742552696,0.00148987960139749,0.00289533462510226,0.0398096824344578,0.0162472241820056,0.00948965076799023,-0.0518638096618763,0.0600174140898688,0.0472304060902588,-0.0355273761393164,0.0202883375903072,0.0234392278378081,0.0422216261339447,-0.0518900550770993,0.0370514884947964,0.0149342752690121,0.0455472842443489,0.0482517623753528,0.0192839286938701,0.05911724067252,-0.00100677106576356,-0.0409349269237998,0.0659956001425141,0.0674249651969527,-0.0747426529204206,0.0431735729416369,0.0108531262403601,0.00681786573268896,4.67374548533606e-05,0.0628869602749529,-0.0176661403186927,0.0786480859660767,0.0524744645035985,-0.0418882346720846,0.00538064945739723,0.0144900607737209,0.0412227917798612,0.0745418645046977,0.0701532715513153,0.000597512655533077,-0.0725604402984308,0.0226924941421917,0.0577535415337205,0.0666303270992964,0.0287497000157599,-0.0541356280325612,0.0614123152109852,-0.0198186871129175,0.0215594750847957,0.0270599850426697,0.0423853336672181,-0.0427899707680521,0.0499074528313328,0.0793047787834378,0.0473748562227119,0.0524015135663986,-0.0379093903107147,0.0793937288558049;
    springCoefb<<-2.43192959884493,6.08514030121006,4.04035678062355,4.0939778826623,0.885759964104128,-4.40632201143234,-3.84532036998163,-5.33669165449445,-0.201902543387136,-5.74421405700816,1.84616840882136,3.98087731087271,1.11369296464456,0.934977743915378,2.02614908115595,0.959073349058476,5.15053093739066,3.006446367785,-1.6732576622727,-1.23768408693475,1.27684178065925,6.04271727318241,3.77882575946868,-6.08085982355576,0.652931759876538,-5.48543869292053,1.81651248891278,5.96869830004477,1.31094790673134,3.23199520554859,0.163817214486977,5.16220350788641,-2.7309598098749,5.48964219698341,5.39277565112383,3.98484668574114,1.44712120543646,5.41086102056708,5.39994485533409,0.946468807561242,4.0832679471894,0.279088417230802,-4.91215907281787,4.6082018532571,0.92838679515286,-1.64664393517429,3.42329263641773,5.70322596816668,-2.65748894621183,1.29342632652648,0.568687125340165,-2.43705856028738,-1.92000280751018,-1.8801132043928,-0.517201690556528,1.68458566290145,4.00263656362903,2.50000791231667,0.452871199207279,-4.12575887624951,1.3802421085177,2.0898833298572,1.09336188408456,4.51445899470899,0.99940722926591,-2.21682713707008,1.5502626395456,5.9937160473333,-3.08915132235467,5.79470203052607,5.86985869733578,1.12879047483547,2.60477802497352,0.454742152576666,0.517470415219057,-5.5002405095043,-0.520598258586578,-5.43747913056799,4.15794818474089,-4.77267751366817,0.0886738094516904,2.80624253330262,6.11985611039204,-5.33693277264827,-6.2580654340777,1.58892567487088,1.83628608125866,-5.35466976461089,2.3785831099318,-3.87135237064131,3.14285901158726,4.07157566778108,4.50171615921588,4.25724858523384,-4.71109470247867,-0.420077554397874,2.01939375860174,3.12235313706693,3.04526415542555,5.21342763624706,-3.62827223889698,-0.737268266315717,1.96095282746985,5.25969098607654,6.00065908626095,3.66813714062273,6.04263567657224,-0.803124378177302,1.71191742881433,0.288838383269688,0.707383314006201,1.04678416576555,-2.49262183233593,3.90141997675589,4.53851627650281,-2.46750206641364,2.16747008949913,-0.845862873766654,2.96766701262683,1.696825034699,1.26885943120368,1.67349275491913,-1.65457166366446,-0.215499085192086,-0.443091929617129,-0.0824811602914577,3.72336076311685,-4.70688337686707,4.86531505313461,1.28618349288566,-5.77664094062,-0.588400267973177,0.716840377031068,4.52819929759633,-1.61189448189663,0.434314257440346,4.15175172698769,2.17904391568207,5.91437508511471,0.300796283741831,2.01864875880472,0.338573193269244,5.38597896424107,5.14448955286199,3.82100555758673,4.50086791031566,3.6169715871486,3.34167081581541,5.25594889042253,3.6242449893384,0.214752047755745,5.92319047857269,4.68127868871957,-3.14173898982101,-2.97035686401395,-0.230860614223697,-0.657414860877198,0.394065895226855,1.19275490160846,5.56181821803563,-2.91342920436725,1.82579686459939,-1.63544369410512,1.79662060059091,3.02465413262939,-1.00109808835509,4.72917326240044,4.64586125535466,2.46719717703016,5.95293752900116,-0.728069876461334,5.10586141337541,2.35273521470591,1.07370842850819,-1.93135481169837,4.37060868490869,-5.19487679371209,0.391702183144271,-1.21653276126522,-6.13711309553413,-3.61258757467371,-4.63071240534989,1.85101678249404,0.308870603002905,-1.4892661951709,5.49822140170814,0.837830493970494,5.16937915605619,5.61545234316673,3.98285231687726,4.89567875026041,-6.23992361757923,-5.34797394185981,-2.9636930323549,1.39715120233631,-2.12644323785466,2.7619378108539,1.80457371889052,4.09084134539069,-4.52002808382579,1.6512736377946,5.63214168372345,3.85473970142402,4.58615899606649,0.42266490637997,-4.35980031027435,0.0279554131341941,1.51097331266788,0.96973753156521,5.09460785186898,1.65704541713375,3.6403351568915,1.05442214408828,0.661927881806177,4.70902626508569,-1.02537075450348,0.192632458299389,-2.51630261939429,1.84969643455691,3.83759113021315,1.5746427306926,5.27170062187801,-4.70822915029589,2.50985398883279,2.6840672558054,-6.16522383252989,-5.89977445487354,3.05056480308843,2.21451166528256,4.47425209051715,5.5492140964325,3.75304063429704,3.82320856838893,3.12076859200485,2.05601442451186,-2.03731173108278,5.0441534817305,-4.19921536235395,5.75684678158511,-0.269294192555964,-5.38779271633664,5.19345331719696,-2.91214424151613,-4.92389726981756,3.36810651670152,-4.48630318228211,0.33391717567896,-4.53575025726799,1.0605535709224,4.47792557723166,-3.00405579583632,3.00532537434112,-0.591793793133674,2.71423009991725,5.14505028459612,-2.34202100170558,-1.31113837281043,5.52846103557426,-4.83914922729897,3.79489135958476,3.12888332563001,-5.57312033671814,1.83308250263569,3.19973161421444,3.83083346113838,-5.96668672335076,3.33493105350654,2.5918017370172,-3.88271688570471,3.38475512344542,4.01203987647953,2.89624094641009,-0.173391873853497,-0.589636697054897,-1.62803646600724,-1.48343160398213,1.64742297071013,4.98906590967172,4.47871139269655,1.19791773758934,-4.04150873496511,3.90776009537591,3.43991137256052,2.61626586932033,1.49626440893562,2.84893529190212,4.7158110969248,-3.53724639953568,4.50491096743946,1.49275078504112,-5.77939771606788,4.17997546079867,4.1568390815018,-6.27213179494185,-5.02828562505671,-3.57958680001224,3.34456559070511,0.127808274665595,0.799359763371303,1.47706369856494,2.56591515626386,3.63946459053563,5.86683490040684,5.73396412750752,4.27641316775727,4.53509332871182,5.61569700322845,2.9393886033297,-2.31855975768557,-5.79668449207325,-5.13523515189644,4.72865749192539,-1.49676813258366,4.1860763159865,3.37062530542857,1.21079109463131,6.04617452801232,-4.80287018786472,-2.17921914776464,-5.35037657410494,5.12690620844648,3.61694724017608,3.05138930367418,-4.01275011252449,-2.72372981137355,-0.0180124155448347,4.12303341309968,-3.74873780023526,-0.887138111904024,-5.41990213502721,-5.9774984601979,0.552717896425435,0.921612258771997,3.17838284590814,1.86797359438216,-1.47430087116196,1.35138812294626,-3.90908982234295,0.382622065522161,3.58663702089682,-1.89477216410266,0.842265657489956,5.32264774649661,4.07279963344888,5.93452553491607,4.58997732369303,0.786062578409914,0.325200986907274,-5.9917002755706,0.588543628002016,-3.55057119553441,1.28408373976525,-1.25341481552174,-2.18636580089857,1.02549154472145,-2.97758908615075,-3.02601767758674,-0.79860588073586,-1.31513570862061,-5.10264541624672,-4.04774808632133,4.46998918241493,1.57982980509186,5.8374249517381,0.290319230405863,1.86919610233567,0.185816115052586,2.65571266392906,3.43118334743,1.54324804109866,-0.19298985557722,-5.47214491078457,2.00854924420448,2.55367598272938,2.67066286720342,-4.40899082454341,-4.04429312445419,5.6046481394415,-1.30713353515645,-3.75280819417312,5.43198984070525,0.316535709131635,3.81446074559214,-2.10461017481648,5.52229986841057,-1.44323290968642,1.20098593903276,-3.78690300917618,4.04134640957772,6.16903543626382,-0.65473740484626,-5.61395256815843,5.11120893505366,5.25665177966897,5.8310232749945,1.77763132207111,0.842662682004633,-0.266345809952913,1.76766768752563,2.41652541542875,4.29474875067796,5.02212468610824,-3.23113956818491,-4.46801644912111,3.88677049527466,-5.81971337358793,2.74959988627578,-3.37151591823711,2.76592180759093,-3.16280122963266,2.80490867323494,2.55004818904424,0.689441449147486,-0.0351284046447834,6.06388826542954,3.55993961157682,0.535486416593954,3.05170925752847,-1.00278238102184,5.46116454437677,5.81344950372511,-1.74189744386769,-0.332778030647898,-4.03931177144349,4.25502161714161,4.25874652459808,1.4159745097112,3.55494616590075,4.0525786915951,5.09901468640174,4.31990406550899,1.04140371424045,5.82363717143333,-2.31027463831456,-0.906303222026946,-0.757257994087879,-6.03794887918608,0.696102823104212,1.00660345861309,4.57369870271486,5.17365517040888,0.38104645055342,3.68186372776027,3.74661889802106,3.58529704079631,4.04925640652263,4.59815072532965,1.95431786572231,-5.33829054200515,1.63610028456171;




        ArrayX3dRowMajor massPosition_i = ArrayX3dRowMajor(massPosition);
        ArrayX3dRowMajor massVelocity_i = ArrayX3dRowMajor::Zero(num_masses, 3);
        ArrayX3dRowMajor massAcceleration_i = ArrayX3dRowMajor::Zero(num_masses, 3);
        massPosition_i.col(2)+=1;

//        ArrayX3dRowMajor position_history = ArrayX3dRowMajor::Zero(num_frames * num_masses, 3);

        double fitness_child = simulate(
                    num_frames,
                    skip_frames,
                    springCoefa,
                    springCoefb,
                    l0,
                    massPosition_i,
                    springtoMass,
                    massVelocity_i,
                    massAcceleration_i,
                    best_position_history);


//    double fitness_child = simulate(
//            num_frames,
//            skip_frames,
//            springCoefa,
//            springCoefb,
//            l0,
//            massPosition,
//            springtoMass,
//            massVelocity,
//            massAcceleration,
//            best_position_history);

////////////////////////////////////////////


    // render this in glfw
    render(best_position_history,
           triangleVertex,
           cubeVertex,
           springtoMass);

    return 0;
}
