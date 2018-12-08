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
const int a = 4;
const int b = 50;
const int c = 4;

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
            glm::vec3(2,4,2),
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


////////////////////////evolve////////////////////////////////////////////////////
    int num_frames = 500;
    int skip_frames = 32;
    int num_evaluations = 2048*64;
    int population_size = 64*4;
    int selection_pressure = 64*2;

//    int num_evaluations = 64*10;
//    int population_size = 16;
//    int selection_pressure = 8;

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







    //////////////////////////////////////////////////////////////// playback ////////////////////////////////////

//    InitialHeight = 0.5;
//
//    int num_frames = 1000;
//    int skip_frames = 32;
//    ArrayX3dRowMajor best_position_history = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);
//    ArrayXdRowMajor springCoefa(num_cubes);
//    ArrayXdRowMajor springCoefb(num_cubes);
//
//    springCoefa<<-0.00807191759258132,0.0057528475186568,0.00892870436372641,0.00816345307510467,0.00326291043890466,-0.00521358444132543,0.00972798142826913,-0.00809914504089353,0.00493460658701817,0.00578044945434771,-0.00653398739943932,0.0087542775546919,0.00521723998487799,-0.00806637220925948,0.00450424979159693,0.00690545989056372,0.00628490383132584,0.00816801807044032,-0.000869716476122717,-0.00131203809348496,0.00430309172925527,0.0068804858452867,0.00155689034283832,-0.0074004802468235,-0.00837469929753556,0.00791695799542608,-0.0072608779451162,0.00846181769634367,-0.00771081937370394,0.00761366348602514,0.00164899256622838,-0.00578273695697204,0.00336651099536871,0.0070624506674149,0.00304167991812884,0.00886444176494351,0.00622778703579864,-0.00755827501302505,-0.00923470328526325,0.000298719066334292,-0.00201400205121096,-0.00576869068470257,-0.000947003388287036,0.00618804367726233,0.00677198934566676,-0.0013248315971926,-0.00989130219439571,0.00299573098914499,0.00925229311829557,-0.000761018661205199,-0.00135096788934943,0.00560501031372929,0.00571864623377037,0.00147829989692247,0.00820453007621902,0.00734394694554803,0.00345952488204518,0.00488298551872512,0.000408271913607808,0.00760937030109613,0.00373950365848752,0.00817217840270213,0.000756988717307994,-0.00770812148121564,0.00189775416250236,0.00996639290205224,-0.00605676687138936,0.00150405502482506,0.00112520021904502,0.00770213844536525,0.00214773347235645,-0.00987756097683569,-0.0090013413964777,0.00379514827594184,0.000436689113470115,0.00109777899975785,0.00978117786337676,-0.00316755265610644,0.00339425745693712,-0.00761552082263656,-0.00278779487721054,0.00847583381388142,0.00845316775071117,-1.4447565197222e-05,0.00348257613996164,0.00145469968740582,-0.00666325863295386,0.00566826890952339,-0.00366231480318229,0.00992899595731506,-0.00111241675499474,0.00311625535279338,-0.00349587821098784,0.00717767871287412,0.00540813387157774,0.00840187595151452,0.00507613014837576,0.00935136700018839,-9.40690329736416e-05,-0.00434330779795689,-0.00399450901616109,-0.00794633556993042,-0.00422086877479258,-0.00299585041263879,0.000521223540660563,0.00621582033867753,0.00810192858711906,0.000302401394724101,0.00315286090278666,0.0089400192997139,-0.00556164922917339,0.000138293081027592,0.00556925467009156,0.00872698388934461,-0.00715762851627433,-0.00410797398263029,0.00122014705614194,0.00289039683197178,0.00780270146490419,-0.00534303764595792,0.0052175929129956,0.00258718631816431,0.00665109326907019,0.00625993746158664,0.00546601892237832,-0.00943094924065794,0.00829447991348418,0.00235164370031638,0.00676105204121277,0.00548864812380106,0.00486090911103749,-0.00846493246893628,0.007600174163282,-0.00654555735017432,-0.00642026714813908,-0.00280428304001888,-0.0011391472868338,-0.00208024142872553,0.00893054976078242,0.00763761644255744,0.00621316308970619,0.00162360799056653,0.000923436270525417,0.00450708348048249,-0.00724252136295779,0.000797576700544527,0.000521238981988859,-0.00464777337137972,0.00475222809927176,-0.00779633661629461,0.00428465982632929,0.00599355980567334,-0.00722632896957282,0.00633334837683167,-0.00626098670822614,0.00268231545564707,-0.000786966975213479,0.00185493257078384,-0.00415136507439956,0.000598360193240624,0.000346795863991597,-0.00621197269121742,0.00341204059003482,-0.00595347204057196,-0.00481309723798795,-0.00677012942580978,0.00727772047645746,-0.00719840145539418,0.00114962914546469,0.00160522608626877,0.00560415060224117,0.000532396691583383,0.00814104028828165,0.000628284125053117,0.000243674185995876,-0.00401986496244551,-0.0044425402276416,-0.00184654428243942,0.001332361675488,-0.00969031213768307,0.000357119101265966,0.00597429388408387,0.0063032476773035,0.00313079013169314,-0.00585486661449767,0.00416473138251347,0.000722291455940478,0.00335816641960208,0.00748337412905272,0.00657092639085414,0.005673551412023,0.00423811955109151,-0.00189465306321841,-0.00613014678290586,0.00675971781684073,0.00660649394439575,0.00689864428376618,0.00258941110809772,0.00728751271449922,0.00374962028061764;
//    springCoefb<<-0.551016871699212,-3.15729872084238,-5.40404032313082,0.242203483966531,1.48689883581391,-2.53950916938092,2.72109785495702,-3.52734512436925,6.02604484579568,3.79031825184668,0.06088809188604,-4.79887660297642,-3.11368640674377,-5.07483487003584,1.90774557220025,2.35465204685259,3.41656503543474,2.54733136108597,-2.00227101667204,2.27033382721358,0.609753700589618,2.56083371829117,4.02997929269357,4.92372114618708,-1.19002872376256,1.87848255837888,4.44436645941365,1.56495141167493,-1.36749622587472,2.62190158028565,2.82472561868147,4.36467210827774,5.74778805944327,2.23028093462838,1.60024210219509,-4.09249610466774,1.24827237455699,4.10604033837472,-1.33665602903699,-5.37596423973794,-2.15394463437049,5.00741726870072,-3.89165564271437,1.55653148026888,0.71748793155196,-5.41951223225552,-2.91297899426166,2.95348895655777,1.53471862090301,1.3679351890663,2.03875023596127,2.26618155667802,1.23190268159479,6.1682775662163,0.906717497013432,-6.24131124801944,5.56013388453743,1.19686486961348,2.72109785495702,1.07917414014591,2.45951810187205,0.998068628430766,3.36554978939063,-3.32338164835877,-1.33107909867942,0.0997814744277313,2.32320655304044,-4.7365890046814,-0.837791515981837,5.68220035847482,5.4147704057438,4.28789925716271,-1.06642538459816,5.42152482649468,2.19376158312559,2.01687455519431,3.13708802807552,4.70881631878979,3.49043729881357,3.07770694366882,2.70532583410795,-1.13323850234129,5.67282346281929,-4.39151004804894,-3.31334182314918,-5.40038648633162,1.63187590029424,-2.86226715639717,2.07966358328187,4.84414020298854,2.40479385577111,5.15797349203374,1.03551146254214,3.09724195893137,-4.4485933621767,5.98761756386272,3.86461511800834,4.15779839671541,-5.03215664667036,-0.285152189162062,4.42411935754224,-5.90057144092656,-2.28043813785103,-2.92549123290759,0.709041636951636,6.19650864527456,5.37456852228672,-2.43705554082451,3.04814519677103,-0.916770131533504,-4.37479512085961,-3.79399983532115,-4.64057809064588,-5.9142018346553,-3.12105698581854,-2.93675051149778,-0.26370585143756,-2.40624018433536,2.64430343516568,-1.82849631179274,2.00196571849143,5.92860423364835,-1.70348644713852,0.26202819456711,3.52093687034608,-0.424073392179973,3.64994480874062,-0.868296371039649,3.70873927887177,-2.83675181158113,2.52957775973524,4.55018831635155,4.12587728483689,-5.4272841248628,-0.650572190411754,-5.67561849064847,5.8309412930667,2.16675671855264,4.24423437749254,0.335600151780653,5.2238411973023,2.48327301288141,4.19371463247202,3.15946310810475,4.91517543429709,4.4255152415988,2.9174478227093,5.8282638077774,5.41431640960629,6.09635058676835,2.53783056453957,-3.55339184783624,5.39686200074684,-1.87781439836787,-5.20240343688496,4.16123025804233,4.37915649505027,2.0828119049317,-6.15118873129875,3.40106295988261,2.67808045098987,-4.46861578909323,-3.16444281094388,4.41490013573475,0.428540459274274,6.19411579117371,3.4195313836648,5.23902014706467,2.07768730387469,-3.50913117430733,3.67941107422957,3.25211395923951,3.16350396031484,2.08407819228096,1.68278765246566,-6.27544037137673,0.862283584728441,-4.32363561251528,5.83600864225234,-0.00658521151694587,-4.5104702315986,3.63623151019987,2.72320814064681,-5.39679343670344,4.95933150918508,3.01146487377684,-1.7870139420138,3.05530279838368,3.19241001079028,-1.65501747331254,2.27851361660629,3.01598643678031,-4.65607616259707,5.91511237240774,-4.01188033752073,2.53480583449695,2.84802611687618,4.32756406829176,1.92334059070792,4.81895540812299;
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
//    render(best_position_history,
//           triangleVertex,
//           cubeVertex,
//           springtoMass);

    return 0;
}
