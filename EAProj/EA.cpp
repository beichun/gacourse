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
const int b = 40;
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
const double CoefaRange = 0.01;
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
    window = glfwCreateWindow( 1024, 768, "Five Cubes", NULL, NULL);
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
    glm::mat4 View       = glm::lookAt(
            //glm::vec3(4,3,-3), // Camera is at (4,3,-3), in World Space
            glm::vec3(2,10,2),
            glm::vec3(0,0,0), // and looks at the origin
            glm::vec3(0,0,1)  // Head is up (set to 0,-1,0 to look upside-down)
    );
    // Model matrix : an identity matrix (model will be at the origin)
    glm::mat4 Model      = glm::mat4(1.0f);
    // Our ModelViewProjection : multiplication of our 3 matrices
    glm::mat4 MVP        = Projection * View * Model; // Remember, matrix multiplication is the other way around

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
//        g_vertex_buffer_data[3*i+1] = massPositionFloat(cubeVertex(triangleVertex(i))+k,1);
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
// Index buffer
//        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);

            // Draw the triangles !
//        glDrawElements(
//                GL_TRIANGLES,      // mode
//                indices.size(),    // count
//                GL_UNSIGNED_INT,   // type
//                (void*)0           // element array buffer offset
//        );


            glDisableVertexAttribArray(0);
            glDisableVertexAttribArray(1);

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



//    int num_frames = 500;
//    int skip_frames = 32;
//    int num_evaluations = 2048*64;
//    int population_size = 64*4;
//    int selection_pressure = 64;

//    int num_evaluations = 64;
//    int population_size = 16;
//    int selection_pressure = 8;
/*
    int num_generation = (num_evaluations- population_size)/selection_pressure+1;

    Eigen::ArrayXd best_fitness_arr(population_size);

    ArrayXXdRowMajor springCoefa_arr = initSpringCoefa(population_size,num_cubes);
    ArrayXXdRowMajor springCoefb_arr = initSpringCoefb(population_size,num_cubes);


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> dis_selection(0,selection_pressure-1);
    std::uniform_int_distribution<> dis_cubes(0,num_cubes-1);
    std::uniform_real_distribution<> dis_coefa(0.0, CoefaRange);
    std::uniform_real_distribution<> dis_coefb(0.0, CoefbRange);

    std::ostringstream ss;
    ss<< dis(rd);
    std::string random_string(ss.str());
    std::string file_name_diversity = "./data/diversity_"+random_string+".txt";
    std::string file_name_fitness = "./data/fitness_"+random_string+".txt";


    int best_solution = 0;
    double best_fitness = -999;
    ArrayX3dRowMajor best_position_history= ArrayX3dRowMajor::Zero(num_frames*num_masses,3);


    //the first generation
#pragma omp parallel for
    for (int i = 0; i < population_size; i++) {
        //std::cout<<i<<std::endl;
        ArrayXdRowMajor springCoefa_i = springCoefa_arr.block(i, 0, 1, num_cubes);
        ArrayXdRowMajor springCoefb_i = springCoefb_arr.block(i, 0, 1, num_cubes);

        ArrayX3dRowMajor massPosition_i = ArrayX3dRowMajor(massPosition);
        ArrayX3dRowMajor massVelocity_i = ArrayX3dRowMajor::Zero(num_masses, 3);
        ArrayX3dRowMajor massAcceleration_i = ArrayX3dRowMajor::Zero(num_masses, 3);

        ArrayX3dRowMajor position_history_i = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);

        best_fitness_arr(i) = simulate(
                num_frames,
                skip_frames,
                springCoefa_i,
                springCoefb_i,
                l0,
                massPosition_i,
                springtoMass,
                massVelocity_i,
                massAcceleration_i,
                position_history_i);
        if(i==0){
            best_position_history =position_history_i;
            best_solution = i;
        }
    }



    int best_i = orderFitness(best_fitness_arr)(0);

    Eigen::IOFormat print_formating(Eigen::FullPrecision,Eigen::DontAlignCols,",",",","","","",";");

    std::cout<<"springCoefa<<"<<springCoefa_arr.row(best_i).transpose().format(print_formating)<<std::endl;
    std::cout<<"springCoefb<<"<<springCoefb_arr.row(best_i).transpose().format(print_formating)<<std::endl<<std::endl;


    best_fitness = best_fitness_arr(best_i);

    int  num_evaluated= population_size;
    savefitness(num_evaluated,best_fitness,file_name_fitness);



//    #pragma omp parallel for
    for (int i_generation = 1; i_generation < num_generation; i_generation++) {
        Eigen::ArrayXi fitness_indices(population_size);

        fitness_indices = orderFitness(best_fitness_arr);
        std::cout << "generation" <<i_generation<< std::endl;

        #pragma omp parallel for
        for (int j = 0; j < selection_pressure; ++j) {
            // crossover
            int parent1 = fitness_indices(j);
            int parent2;
            do {
                parent2 = fitness_indices(dis_selection(gen));
            } while (parent1 == parent2);
            int st, ed;
            do {
                st = dis_cubes(gen);
                ed = dis_cubes(gen);
            } while (st > ed);

            ArrayXdRowMajor springCoefa_i = springCoefa_arr.block(parent1, 0, 1, num_cubes);
            ArrayXdRowMajor springCoefb_i = springCoefb_arr.block(parent1, 0, 1, num_cubes);

            springCoefa_i.block(0, st, 1, ed - st + 1) = springCoefa_arr.block(parent2, st, 1, ed - st + 1);
            springCoefb_i.block(0, st, 1, ed - st + 1) = springCoefb_arr.block(parent2, st, 1, ed - st + 1);


            // mutation
            int mutation_point = dis_cubes(gen);

            springCoefa_i(mutation_point) = dis_coefa(gen);
            springCoefb_i(mutation_point) = dis_coefb(gen);


            ArrayX3dRowMajor massPosition_i = ArrayX3dRowMajor(massPosition);
            ArrayX3dRowMajor massVelocity_i = ArrayX3dRowMajor::Zero(num_masses, 3);
            ArrayX3dRowMajor massAcceleration_i = ArrayX3dRowMajor::Zero(num_masses, 3);


            ArrayX3dRowMajor position_history_i = ArrayX3dRowMajor::Zero(num_frames * num_masses, 3);

            double fitness_child = simulate(
                    num_frames,
                    skip_frames,
                    springCoefa_i,
                    springCoefb_i,
                    l0,
                    massPosition_i,
                    springtoMass,
                    massVelocity_i,
                    massAcceleration_i,
                    position_history_i);


            if (fitness_child > best_fitness_arr(parent1)) {
                //replace this parent
                springCoefa_arr.block(parent1, 0, 1, num_cubes) = springCoefa_i;
                springCoefb_arr.block(parent1, 0, 1, num_cubes) = springCoefb_i;
                best_fitness_arr(parent1) = fitness_child;

                if (fitness_child > best_fitness) {
                    best_i = parent1;
                    best_position_history = position_history_i;
                    best_fitness = fitness_child;
                    std::cout<<"best = "<<best_fitness<<std::endl;
                    std::cout<<"springCoefa<<"<<springCoefa_arr.row(best_i).transpose().format(print_formating)<<std::endl;
                    std::cout<<"springCoefb<<"<<springCoefb_arr.row(best_i).transpose().format(print_formating)<<std::endl<<std::endl;
                }
            } else {
                int parent_for_replacement = fitness_indices(selection_pressure + j);
                springCoefa_arr.block(parent_for_replacement, 0, 1, num_cubes) = springCoefa_i;
                springCoefb_arr.block(parent_for_replacement, 0, 1, num_cubes) = springCoefb_i;
                best_fitness_arr(parent_for_replacement) = fitness_child;
                // put it somewhere else
            }
        }
        num_evaluated+=selection_pressure;

        savefitness(num_evaluated,best_fitness,file_name_fitness);
        //savefitnesshistory(i_generation,best_fitness_arr);
    }
    //std::cout<<best_position_history<<std::endl;


//    best_i = orderFitness(best_fitness_arr)(0);


    std::cout<<"springCoefa<<"<<springCoefa_arr.row(best_i).transpose().format(print_formating)<<std::endl;
    std::cout<<"springCoefb<<"<<springCoefb_arr.row(best_i).transpose().format(print_formating)<<std::endl<<std::endl;


*/




    // playback

    InitialHeight = 0.5;

    int num_frames = 1000;
    int skip_frames = 32;
    ArrayX3dRowMajor best_position_history = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);
    ArrayXdRowMajor springCoefa(num_cubes);
    ArrayXdRowMajor springCoefb(num_cubes);

    springCoefa<<0.00611268440080466,0.00896622392250863,0.00595924631861012,0.00910507312911025,-0.00954917111413049,0.00715481074394417,0.00852123820742361,0.00880657828943631,0.00923747488393857,0.00658075965843751,-0.00693954593359472,0.00771913961401169,0.00554680250625692,0.0017962372616392,0.00720807363475542,0.00882129560169824,0.0080975856271183,-0.00960494529437504,0.0048715992899945,-0.00845857379886255,-0.00167299557089479,-0.00972984191949006,0.0074788131725541,0.00562413189438021,0.00681799093620638,0.0098819015125528,0.00621981708335722,-0.00975585160765604,0.00997701433317959,-0.00844637337068393,0.00641493455663543,0.0094981806071001,0.00890318147859004,0.00899943768435204,-0.00955879997441489,0.00838131389880351,0.00945687374270121,0.00803922885730423,0.00892535316707816,0.00924448879527596,0.00498607140298981,0.00707545620253098,-0.00495422852921962,0.00472021113835285,0.00610129088913151,0.00802675297305993,0.00987785590493843,-0.00990346824745809,0.00872619908243706,-0.00730446984865911,0.00838895282054516,0.00268168829040681,0.00881743827941107,0.00896024449869608,0.00814350408415473,0.00758983390854558,0.00818012204292337,0.00868827508199779,0.00872845261817176,0.00629767122857552,0.00840018532697329,-0.00756185261419129,-0.00988898718724446,-0.00670536066252057,-0.00633483850226497,0.00866399886645502,0.00307570119065964,0.0088723096348536,0.00639012198789095,-0.00770094042536846,0.00894657976040144,0.00812113110368831,0.00655227257311187,0.00244164233597392,-0.0044460443567699,-0.00895347068037534,0.00431573790233384,0.00544633091077158,-0.0088569389185202,0.00590910891702999,-0.00613732353604274,0.00724041402463504,0.0088467648771752,0.00355458150132021,0.00877330077901344,-0.00649026886396588,0.00920440858891398,0.00276702802291467,0.00914405084933351,-0.00990875200364215,0.0067518912630649,0.008578903305821,0.0099238816555235,0.00461397996852825,0.00977396192925268,-0.00341070631212122,0.00290583411832612,0.00864792830247801,0.00994045504980913,0.00812395446021697,0.00846754038396279,0.00799734064365665,0.00873485927201258,-0.00635477876120935,0.00307079424366882,0.00926017795204502,-0.00475244983786365,0.00518841762028986,0.00726970787153158,0.00972726448077442,0.00544213638452399,-0.00542492702390297,0.00630396715047553,0.00898135861522617,0.00920745529847777,0.00908834244744227,0.00501248023240477,-0.00531672338737022,0.00647855006450362,0.00196409406531271,-0.00763194424921272,0.00219306528204729,-0.00682653594614311,0.00479775577254864,0.00917929733851443,0.00615043702817951,0.00681313998197743,0.00897162292173789,-0.00603893536889876,0.00272758377329084,-0.00238281828927939,0.00897281931665392,0.00721207766663845,0.00953679756863435,0.00636170582022578,0.003289720695093,0.0084065170300579,-0.00655284104708249,0.00804457621353092,0.00788136462582339,0.00340355066906612,0.00723698508133249,0.00544255035734287,0.00312984340204424,0.00967360390535739,0.0036120044552268,0.00753244429711832,0.00719906789270427,0.00522440345735494,0.00471331715337621,0.00733393045982878,0.00900754919229427,0.00484967048071671,0.00785956237035276,0.00707877229763138,0.00964719540423118,0.00219721764198894,-0.00123041973972247,0.00655221637646828,0.00245338976630758,0.00399266047411995,0.00236498475650557,0.00687437106732278,0.00843141356415249,-0.00992347133342338,0.00827326623642504,0.00648795390422422,0.00329456366415379,0.00388215481425579,0.00812330109060397,0.00476208593452446,0.00370594073423929,0.00685099295666407,0.00739267313731493,0.00971345753976165,0.00750619794865724,0.00692269462017468,0.0081719201606808,0.00838341682018677,0.00480039815176297,0.00818032519807123,0.00839376378522951,0.0074549775422941,0.0050446057017169,-0.00471497433945303,-0.00650864571170353,0.00659488163122561,0.00796987666068975,0.00555138078116681,0.00322032923959722,0.00658217503386455,-0.00423830488428395,-0.00419745495272682,0.00277022448097381,0.00909030233467478,0.00998330317620094,0.00850767525774784,0.00871666779681885,0.00914089606768512,0.00495785015452335,0.00761020216048239,0.00788411448611138,0.00813628674612132,0.00860746907503518,-0.00767389707158967,-0.00695159351311233,-0.00267251588994289,0.00405487390566067,0.00725932856894067,0.00921913442165551,0.00289058428857968,0.00985180724710433,0.00237231445289545,0.00592657480061881,-0.00989461738611321,0.00668535920636978,0.00938825936007143,0.00710289272877212,-0.00887923685781622,0.00578387118771207,-0.00868908693952909,0.00550486171330917,-0.00140582506144691,0.00476944178238604,0.00727407618406812,0.00237220391709151,-0.00692288777647674,0.0030273709049099,0.00560137980625601,-0.00602418056504064,0.00627060530533577,0.0056078286113875,0.00460098257428879,0.00864507956117154,4.4827391041728e-05,0.00998424071539665,-0.000906628999349954,-0.00536746551998307,-0.00560806127991903,0.00786660582856278,0.00339391790760391,0.0037981592855289,0.00321468436686229,-0.00958411993439501,-0.00668752183983453,0.00533926576158929,0.00369081179656655,0.00336898808990092,0.00438317251486467,0.00822200241415855,0.00705362630612643,0.00850631750871978,0.00847085405336268,-0.00520435715802217,0.0072719419700831,0.00437977540272814,-0.00417614512805647,0.00535131050137985,6.74845290391201e-05,0.0051227979432432,0.00497965788234941,-0.00795925434583763,-0.00796805683428797,0.00688358145738685,0.00336989147282603,-0.00808549564242619,0.00697664081943319,0.00978916643150472,0.00966117817064513,0.00885778237510452,-0.00318773491922195,0.00711209801092735,0.00571840644431045,0.00502333872024084,0.00892191338017691,0.00957483756567037,0.00828771074411667,-0.0089097224170853,-0.00415268016706811,0.00989050558965833,0.00746909848188815,0.00289685175386908,-0.00930807542954948,0.00325322098271957,0.00910940679190931,-0.00223496650915359,0.0082757968850764,0.00802947782819601,0.00941251441898824,0.00632416649116025,-0.00319438706766553,0.0049994939060775,0.00839458791878538,0.00990905971033594,0.00969886917042742,0.00306447503513446,0.00571516283793142,-0.00208175605725579,0.0022418828365588,0.00156306800970951,-0.00797912732603919,0.00974128935134288,0.00440626329257474,0.000526963897705194,-0.00762740324140871,0.00711114906435178,0.0080153868058768,-0.00888661608979414,0.00970681938511067,-0.00613729337050454,0.00914857408923496,-0.00846216184946809,-0.00840456181131516,0.00984049865968549,0.0087923764334952,-0.00911814509849909,-0.00514522040968073,0.00662919222919711,0.00876058861892348,-0.00816390073772701,0.00940895287494079,0.00256418048039576,0.00948816185692356,0.00943745223779113,0.00586530582634781,0.0082733314714988,0.00918505128385361,0.00233352650531266,0.00272513084839396,0.00663868822395203,0.00797300698231091,0.00924257524259957,0.00948283297917006,0.00424972892471111,0.00880135615874761,0.00462778373116201,0.00994854042304146,-0.00951970269881175,0.00218006389006492,0.00543079092060907,0.00967279913563162,0.00685533704865524,0.00887109935622966,0.00823198905493733,0.00771069109333246,-0.00557950492742448,0.00819031455202137,0.00618393333451074,-0.00981753271995929,0.00953943730288262,0.00802057286632274,-0.00601817910374058,0.00957784884589624,-0.00844931540938528,0.00330174239971846,0.00565139586685923,0.00181837640321738,0.00287701290899979,0.00380812233124063,0.00948972859286259;
    springCoefb<<-6.13707627096085,1.44306050185137,1.716604031868,6.01067536441671,4.2052448844477,-5.98337528875675,6.2707392509959,0.844216502032975,0.559353106511185,6.08409864084948,-4.64637490929046,-2.51968455239212,6.04635741597282,0.236723754948917,2.16472408797332,1.09382157389914,0.229615225640349,4.03333459840942,-1.58532090897299,2.22387611409131,-3.58325552931724,-4.90308031455035,4.6990768239357,5.89716952451594,4.22322340326946,5.05762056722683,5.31927877459654,0.552077756466005,3.98645794884991,1.24015868889651,4.22821961127945,-3.04315393687926,3.59304676679195,3.78328623850418,5.6105675805877,4.04499120621868,2.46000661527463,3.00652348961774,-3.85868727461365,2.59959041499554,1.93419553402411,2.75399248912486,-1.19184950668384,2.20927891878173,1.41149070683435,1.1027107342901,1.65119106393904,-1.45141198505346,-5.70884059895999,4.4090636895659,1.54807489946366,0.465241784634279,0.454664758420383,0.514555352798606,-0.432817744935137,0.611051080719834,0.624828789745492,6.26735935053486,6.21540148736463,0.879386165272471,6.12580901270167,1.82829702139194,-2.70258683604027,2.03776711991828,2.12326491379173,5.34059250362993,-0.769983888989908,5.18182034126483,5.07005754048825,1.16296466367344,4.89025809794199,4.36000862209611,3.69954811036378,4.87579011240952,1.30253907130625,-5.13139846204371,3.83413247594085,4.00662012069817,-0.29962524709717,3.61286458738738,0.262005653925718,2.4537259887019,3.21733997514327,3.96893034740153,2.40900348631748,6.01119168663985,2.63970323746843,2.40527631444231,2.38265876464859,5.76502599739228,1.88948201690375,1.63701329626333,0.984170610875712,-5.02719679270339,1.93597468186022,4.56485631689499,1.2300136330654,-5.43268185946199,0.592434242444349,0.518492297759169,0.596769916430934,5.93424957502971,6.24290559677194,-3.4158927240529,4.79781133841148,0.0568492410519931,-4.32450711961402,5.5111231481525,5.76825343809685,4.77892008835097,6.10065002202838,-5.36147564561309,3.073094723745,4.48851636555909,4.87758054964818,3.99887693473515,4.46575023105878,-0.17030261224175,4.54223338891581,2.41935017524292,-4.02790112179632,1.48536805105651,-0.793411245976088,0.953171366701376,0.804763145515785,1.25677454326295,0.464475098056513,6.08144194341415,-3.25620166002581,4.64206167215768,3.71323847271763,-3.82004227536408,4.82687890263428,2.18976315610684,3.12450501068769,5.94354389727728,6.17799260041048,4.23043962199448,4.8243337625403,6.06226046480641,2.15457777848791,5.40029805691714,6.11426548098015,4.65623268823929,4.64613280828748,4.76092215347378,4.37103788501228,3.9678119241622,-1.73671739732676,-3.39939970397749,2.6782651532786,2.95511781048191,4.01881520303059,3.44264459194946,3.85622540589064,-1.95739471211555,1.32330681812954,3.97039161867743,3.16583921223597,0.937081291517415,1.9665385742042,4.34693331985175,2.48642878218639,2.18494528273141,5.99118701954851,-4.56069832927568,4.29074252436642,1.5001183993084,1.52729683185275,1.55676868007094,1.14062983729742,5.61012705773358,1.86054870103106,-2.74858212076536,0.750603604435709,1.24397462512567,-5.48423310800308,0.117235246668236,4.79632711990699,-0.22945188265346,-5.7515403038142,6.26245504022549,6.22495475099375,4.6017691685024,-3.27447054558654,3.12602933895736,4.56470486044463,5.65343733180678,5.23373091249888,4.20255966368648,5.86748491424973,-2.54868468504859,1.7678173147843,5.65423723024666,4.74546057279063,4.52556021925102,-2.63730058342616,4.66762434994811,4.30333007250761,4.10781901008368,3.83961225683623,4.04208478528178,3.24724223068869,3.37080661506469,-1.49083013064342,-5.2631804259964,0.111287520674044,2.4026150440541,1.64426133320317,2.5856362747261,-3.91888405286946,2.29379471585706,3.20338131533214,0.809022738547171,4.16075240862714,1.48911503477411,1.35351866899091,1.50329770904543,-2.4790978103048,2.40144976266797,-2.22657595276849,2.4362311223323,4.91125079966257,0.508404524765036,1.53922244236058,1.2452337257651,-2.10727284367172,3.0675269148287,3.68924876888052,-0.176534514419235,0.248958829912458,4.44299366830365,3.81294911170103,3.4012779288367,5.64090092705198,3.12157590264213,0.839016411417611,4.24290095694377,5.09199699371957,1.07804367035448,5.78857172789337,2.68236887347714,4.25120330670707,-4.54099571380099,5.00749039704845,-3.13891595588858,4.1811767995416,-0.118471226596468,2.12356945548597,-5.53097848933168,4.08071651231521,3.02959689679131,1.13630624141227,5.96952985192559,2.95213226440463,1.68363799442845,-2.73477467884193,1.58713173846586,0.334008272998607,2.40038260192645,2.17003439230653,6.11014738450001,5.09196693953105,2.93234183980809,1.69433772928714,-0.946560772441851,0.595442237970135,1.46475240185139,5.72636602006221,1.14765458925757,-1.89941708692241,0.712920466486293,5.13456878091712,4.74435612609501,5.42207531972464,4.78800202705352,4.87558768972228,-3.43670411437983,2.67821505426782,4.98313762557938,5.25865712736986,5.62923039176899,0.414567110740063,5.47090974760878,4.3004234853171,-4.54340731665963,5.49060152610847,-2.63415497606261,4.48060265264943,4.17874845018313,2.35628686363072,2.85656170522079,4.14527459773034,3.64449689245516,2.98674230126173,3.61790824674415,2.92584608469719,-2.4342772663329,1.88392938998607,-2.5409476976661,0.778481546947622,1.36627659714415,2.07940773391646,2.67465407353201,4.75973611338162,0.884753539743558,1.17484897107679,-3.25333860650707,0.391032685899513,-2.43012118050705,-6.03552523358601,4.22475897659848,-2.27818821669582,0.66222707715405,6.238534238932,3.20046513751784,-2.77673454507402,5.73246843517744,0.0397998318895946,-3.97908175568948,1.02082631646848,5.97822723834167,6.03147384002345,5.88786194241009,6.22597842139754,5.67438824641577,5.86556635714366,4.41738121549766,5.72407691009153,4.14820723318761,3.45966230283058,4.62904215473217,4.1697790870119,3.42525017406157,3.06474821295591,2.34693051501917,2.18094281633376,5.46400687904648,2.48243108370983,2.62869050578458,2.18925631802916,2.29223712616404,2.75849930176571,3.11177327559799,3.19635958673329,-0.625611303241902,2.60471309634644,-3.20112526466097,4.52764866974448,3.0189951505435,-3.50870256857557,4.96368717799533,1.75069706566697,-1.6011561542775,1.42449799955208,1.17039142375124,-3.32930504452514,1.64908921810256,2.88616640462439,1.96055057424255;

    double fitness_child = simulate(
            num_frames,
            skip_frames,
            springCoefa,
            springCoefb,
            l0,
            massPosition,
            springtoMass,
            massVelocity,
            massAcceleration,
            best_position_history);

////


    // render this in glfw
    render(best_position_history,
           triangleVertex,
           cubeVertex,
           springtoMass);

    return 0;
}
