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
            glm::vec3(8,10,2),
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


////////////////////// evolve ////////////////////////////////////////////////////
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
//
//    InitialHeight = 0.5;
//
//    int num_frames = 1000;
//    int skip_frames = 32;
//    ArrayX3dRowMajor best_position_history = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);
//    ArrayXdRowMajor springCoefa(num_cubes);
//    ArrayXdRowMajor springCoefb(num_cubes);
//
//
//
//    99350353935,0.0699945347709556,0.0279196183886005,0.0372754298696646,0.023275334808638,-0.0353948253185464,0.0667176000712056,-0.000708868392141939,0.0436922640929428,0.0607062414198677,-0.061920832312629,-0.0489750701789628,-0.0255419144525854,-0.0272161411993281,-0.0758359912763517,-0.0636525478137902,-0.0389819724666802,-0.0668556268451994,-0.0471462438847619,-0.0513265106693499,-0.00158488024099957,-0.0505461882662709,0.0335339989296785,0.0455900448400481,0.0689644761332145,-0.0487980246398589,0.0454743329274814,-0.0760199969988409,-0.00590198915726598,0.0539418101934445,-0.0422212101343187,-0.0232120541963782,0.0439363448898943,0.0656984083287876,-0.0659366244012195,-0.0127883203014677,-0.0496964170642646,-0.0792190243300139,0.0665028113063904,-0.0225080636621025,-0.0414267164289144,0.0412245066283385,0.0663494203176114,-0.00841421615724182,0.0420215482832964,-0.0233024226610094,0.0130916214422749,-7.7103413677353e-05,-0.0743279516670517,-0.069005969236142,0.049077304140235,0.00880122866891382,0.0559582856930598,-0.0653386633029853,0.0633224085640732,-0.0575099508359609,-0.0174339168227436,-0.0646385323184722,0.0388923018979338,0.000972825177466881,0.022208964201719,-0.0452487389581505,-0.0649152328003735,-0.0619174795513588,0.00733346116139248,-0.0328625364382111,0.0565694240371554,-0.0141776890746214,0.0687525449640828,-0.0545089425027878,0.0433142473377819,-0.0526741715393375,0.0667155642000565,0.0296636675808875,0.0189116123034207,0.0287371124088471,-0.0736387550801219,-0.0479967662543043,-0.051339991079336,-0.0679667067471737,-0.0370027354904463,0.0777373131354048,0.0208345219217401,-0.0610444498718923,-0.0676013502420864,0.00415693041130757,-0.0385544007078532,-0.00503526706482995,0.0195183980928354,-0.0796620988844252,0.0759375581871427,-0.0382726377054456,-0.0449108378425757,-0.0689776746877365,-0.0201901172568044,0.0424226233188168,-0.0218402111259476,-0.043620693219649,-0.0517550658303104,-0.0330876662363706,-0.0181296357224368,0.0715591815074716,-0.0057618377757081,-0.0314140715968861,0.0212228490883591,-0.0668502255467932,0.0773230408864669,0.0275840940082372,-0.0348469918010975,-0.054016950267375,0.0396173872610635,0.00815027270845617,-0.056279637206476,-0.0195480908917021,0.0271058228365638,-0.0438809873740566,0.0646088395941112,0.0685514221287106,0.0310837455611135,0.00412723761244083,0.0688893232442855,0.0270213036737504,0.0458545999069952,-0.056021514672796,0.0380436289860139,-0.0543355173498092,0.0664011086460208,-0.0637965822144396,-0.0179562105694582,-0.0653539571842896,-0.0168842484508102,0.043914153708105,-0.0737947757513238,0.0573539137734817,-0.067499917888781,0.0274280733370353,0.0705036883011943,-0.07017687707682,-0.0249878327292334,-0.044343303574409,-0.0441938273441949,-0.0653704455426757,0.0438069692085529,-0.0204734644761651,-0.00491853635987199,-0.00908720802938901,0.0156455481497783,-0.0203096968402666,-0.0205357859006784,-0.033270706363614,0.0638175407721743,-0.0316464627308987,0.0737505973846422,0.0296721406791695,-0.00766797732918895,0.0317942262961502,0.0553366233293603,-0.0212668687576739,0.0479976441562165,-0.0426195873146037,-0.00662082594196351,-0.0488866043690995,-0.0787054336810044,-0.000415601618781496,-0.0715326905956178,-0.0662053515697854,-0.052987528356252,0.0789709977055765,-0.0563822286466054,0.00202463898902044,-0.0453723059433383,-0.0205760559162945,0.0166541934463448,0.0784346632652146,0.0389504796075404,-0.068264342988033,-0.0106525447641744,-0.0254039723171871,-0.00857403982829957,0.0488116693351472,0.0213253213191989,-0.0247564990561253,-0.0628347934702573,0.0150759186293352,-0.0750843584514616,0.00949722920055372,-0.0331298550745146,0.0602522648778987,0.0682303605173856,-0.0651322109928039,-0.0623673225112107,-0.0183904654990837,-0.0340188153619034,-0.0610727561922151,0.0611939328821348,-0.0255515059575213,-0.0472781076874948,-0.0717935954741172,-0.0265805083264506,-0.0236603363341002,0.0102310435149032,0.00804718580471686,0.0357636077496054,-0.0531147631132578,0.00648184899542567,-0.00528591271736003,-0.041379106026785,0.0758293042312513,0.0493101149654529,0.0300468541449154,0.0446409734918927,-0.00936456371534829,-0.0747096449857157,0.0618061800961411,-0.0742886451605189,-0.0697940034371773,-0.00869659077781094,-0.0274185002350335,0.0704582614407215,-0.0204662303349312,-0.0125507112278374,-0.0719090611449951,0.0411433042404909,0.033430473484765,-0.0529818172627044,0.0223372370481199,-0.072121032547262,-0.0202599249501992,0.0305436415740026,-0.0187015407992068,0.0360797387157007,-0.0392253149856,0.06934564500551,-0.00815665353469395,-0.012340078024352,-0.00417250607357011,0.066557433747946,0.026280815948863,-0.00834320184231885,0.0358675487133989,-0.0236723299807274,-0.043702228424932,-0.0534970150764552,-0.0183819748919373,-0.0618960484032966,-0.0477856602369741,-0.00817597832911461,0.00940736089339822,0.00479583960249826,-0.0177177168883931,0.0689411305584671,0.0722451283746609,-0.00962677803338821,0.0300844347244522,0.0256756017849201,0.0173914047039074,-0.0275783282274279,0.0335545693121639,0.0771314797537082,-0.0770346867279311,-0.0651469715615488,0.0332112184694089,-0.0362600016390253,-0.0758013266305445,-0.0549454351397908,0.0313999203366227,2.61672958853509e-05,-0.0683880014663506,-0.0223192637890201,0.0716829654535665,0.0474795472470483,0.0340084062302524,-0.0520192629713655,0.0739825321705931,-0.0643735687361907,-0.0339153113746621,-0.053803128066381,0.00745045300920049,0.0554920495187361,0.0309927115361172,0.0697327361208074,0.0444331800772032,0.0232378398362723,-0.0198940419870867,-0.00548238527285046,-0.0310865583788075,0.0774973627913265,0.0469392864997216,-0.0775319891411494,0.0746288424705289,0.0499045998462963,-0.0626789607026982,0.027840060865432,-0.0663554018672348,-0.0584802873332427,0.0528946257256412,0.0450445184693879,0.0215458799626426,0.0645066243337964,-0.057274745394138,0.0132288454162091,0.0319861715063388,0.0567336608361144,0.0412095824448436,0.0259687036769319,0.0723600921744295,-0.0727057290043243,0.0521655756105509,-0.000189454890875824,0.0627863205889176,0.00315828707216228,-0.0104567188445743,0.027219500591615,-0.0536038731660712,0.0496492391683391,-0.0582628847557413,-0.00469043154487872,0.0471466018851598,0.0686764018184861,-0.00222242068602816,0.0417754443556887,0.0385810015902766,0.0150986186112736,-0.010384494853385,0.0522255997230418,0.0366183312780309,-0.0374898691277438,0.0172701181179239,-0.0218357887593265,-0.0529832448684532,0.0399953727237859,0.0713930566568827,0.0590029266378856,0.0167290334853945,0.0326026390272205,0.00497163024031167,0.00908912558531814,0.0398969100228962,-0.0228627942236433,-0.0711003293800635,0.022683230537308,0.060295492848519,-0.00155704822463777,-0.0300972689455828,-0.0733083803175522,-0.0319078090562987,-0.00836015362681828,0.00200118813756909,-0.0647612072456447,-0.019683751882838,0.0797787675260467,0.057014237110044,-0.0611027503670672,0.0148773860628146,-0.033370257743341,0.0711228493559746,-0.0285042826591545,0.00913987312891514,0.00839296739939274,0.029659928581519,0.0361566282604619,-0.0316116599513272,0.0210529851638959,0.0151595548238417,0.0651173735340672,-0.0263443758833894,-0.0598688149358466,-0.00579350088061462,-0.066447465934999,-0.00273160915948991,0.00310616973932188,0.0362357646768148,-0.0224361163109709,-0.0784508784853159,-0.0738615043432738,-0.0157444966285231,-0.0303586875416146,-0.0022216579700921,0.066256691509046,-0.0151198947872594,0.0580945901470699,0.0660354589605869,-0.0381056576772154,0.0769918397800028,0.000912845023401481,0.0085240845794436,0.0681146890614716,0.052408562364247,-0.0623360423661471,-0.00349234361364149,0.00206849087126017,0.0538205858943148,0.0448959964350313,-0.0568785240393498,-0.0110198593563493,0.0300133699690985,-0.0032228999227392,0.0091113257078041,-0.0557801309860219,0.0103296341422618,-0.0736202834516858,0.0273260387533,-0.0334346012554293,-0.0160563997626567,0.0288751602679841,-0.0272961055241973,0.0481991036088202,0.0775463980135299,0.0504822365057106,0.0344557950433603,-0.0166034220608898,0.0285768266527806,0.0204912540039473,0.0252909202618948,0.0255686663582775,-0.0585959010471571,-0.0461849952331674,0.0136833554197491,0.0738126613170899,-0.0285210375993145,-0.0698089882683982,-0.00411884788615575,-0.0547004517795055,0.0550870082411389,0.0190026280744945,0.0142796888641453,0.00510037813573161,-0.0642202719227505,-0.0566089854280506,0.0293202471497097,0.0261093622195112,-0.0502292688797364,-0.0233537141714961,0.0726747610385878,0.0137143313576068,-0.074478553903512,-0.0346213445601153,-0.0180865651080788,-0.0759620812516483,-0.0641391081289105,-0.0636307701392243,-0.0125655033125381,0.0444377185238701,0.036860483864723,-0.0672745831251492,-0.0099936151178524,0.0582645828175659,-0.0334595783583166,-0.0763102597726091,0.05207724406015,0.0180193840423689,-0.0661192480410073,-0.0320416039005116,0.0433189322628635,0.0689677602001316,0.0669610241739829,-0.0224013788729912,-0.00593186173864261,-0.0772592478232734,0.000989635698958136,-0.0566116146634387,0.0288501144707436,0.0307603668192217,3.46712395710291e-05,0.0215248754348256,-0.0355253018976773,0.00555611733605904,0.0669035308747103,0.0263881329942439,0.00959403608441075,-0.077235577328706,0.0427573628550197,0.0770285328463784,0.0472021412696699,-0.00038215328025732,-0.0702460503532766,-0.0427914739226883,-0.0221175705371972,-0.0237056286370874,-0.0391017336952974,-0.0500403264770472,0.0743137554052816,-0.0252209817363047,-0.00208193037755877,0.037632687668145,-0.0362532215361731,-0.0151209062035759,-0.064768691279352,0.0378149167251843,-0.0123801540268493,0.0162209444196061,0.0612033021362514,-0.0635300396306114,-0.033018688835678,-0.0187620266986834,0.0379948358042142,0.0114560092666447,0.0667940906373756,0.0248983666044187,-0.0421558578136172,-0.00361187335271941,0.0276627893502185,-0.0793984950331033,-0.00658334058084681,-0.00513506945461737,0.000219351686639406,0.00317060914038244,0.0320734566226944,0.058101781223948,0.0594649805032951,0.072971722927397,-0.071938545327605,0.0537787358340708,-0.0322492588834135,0.00597952429483623,0.01141142342771,0.0114975195804134,0.0708586180912604,0.026642732148358,-0.0306875636944024,-0.0215215359355889,-0.0371363235065417,-0.0494842616326568,-0.00505157556620034,0.00984498765778028,0.0117537116686598,-0.0470567397619862,-0.0586990031500808,-0.00145219776847036,0.0578416268424325,-0.0208548609636979,0.0749359288788102,0.00550441611814518,-0.0202533559222954,-0.0116474117020366,-0.079630653410978,0.059965995764344,0.0715231974383459,0.0324428032117164,0.0380677769137862,0.0509881778671351,0.0254145260646076,0.0461292315861812,0.0247669136267001,0.0731652672556999,-0.0278912441189826,-0.0438216629455898,0.00466278676160741,-0.037032626102228,0.0628210692027682,0.053975223067205,0.0214458379621831,-0.0543152543037735,-0.0755090386399576,-0.0636057376040172,0.0355297333540068,0.0162446730287022,-0.0306624773660034,0.0568307302784318,-0.0652075248142739,-0.0528208505980768,-0.044024130759772,-0.0702715959354637,0.0326835655945742,0.0157225133179326,-0.0019190076375003,0.0330529121835962,-0.00431149091772339,-0.0103958102736603,-0.0145042846046874,-0.046243714078443,-0.0394076324810309,-0.0690897586145856,0.079885517582244,0.0653592811456692,-0.0759244914333916,-0.0280057266112444,-0.0584623818744264,0.00873829540272164,0.0149616472865276,-0.0756413126716583,-0.0172864816045791,-0.0435925147512893,-0.0499565669754318,-0.0127955201700309,-0.0271982523553065,0.065573166378575,-0.0765508472158345,0.02213927027869,0.042403896582501,-0.0617583719556026,0.0493184196806133,0.0783797658972348,-0.0520299678910663,0.00200198520068172,0.0141022791406616,0.0260510244714334,-0.044945102615722,-0.0702092118515676,-0.0643447858767327,0.0205506127795906,-0.0364529258555048,-0.0237524183577636,0.031460854165005,-0.0365674083477666,-0.0383931372866003,0.0355363628061192,0.0154268650409891,-0.0168555190865209,-0.035725341865665,-0.0496114877469891,-0.0124968317581792,0.0269881765297559,-0.0132040024982784,0.0175466012663891,-0.0658073437147808,0.0395977451464151,0.00311976757045825,-0.0623581908561095,-0.0182629846494007,-0.0344763358470408,-0.0441165628117121,-0.0489445649687874,-0.0360965700243118,-0.0161465307027784,0.0330574202318943,0.0580057091163498,-0.0700955063058508,0.0681123176161723,0.067796497339288,-0.0544402921825835,0.00866293032125707,-0.0486564285907226,0.00180728945965287,-0.039876215513738,-0.00522383693848916,0.0434141521730526,0.0756601472923812,-0.0697969719720059,-0.0534413669879741,-0.0400651946477896,-0.039408459718995,0.0140618012538468,0.0669229818819663,0.0273875378572324,-0.04839159755427,-0.0788843618328145,-0.0130147170708583,0.0347281700906941,-0.0612425526889239,0.0487222983542468,-0.0797481658308525,-0.025359115500636,0.0797777333854594,-0.0358447358551643,0.0384943537965856,0.0328351535428479,-0.0578390268133204,0.0483988474907348,0.0209474711590202,-0.0700425295485382,0.0739585553081513,-0.0503895985942285,-0.0386989581392607,-0.00423415530670162,-0.0102658141079665,0.0360772049222501,-0.0408200032081548,-0.0146056668900911,0.0462802329502442,-0.0142613701216231,0.0253291384621193,-0.0731282267687508,0.0798004311322237,0.0122521203440857,0.0342593110884816,-0.0485911664220463,0.0133677585112712,-0.0587554060568825,0.0661370036686477,0.0321252058223473,0.0699668922973643,0.0663888378377952,-0.0732339097527945,0.0697446256083179,-0.049455898091875,0.0452604440437911,0.0225797791511657,-0.0272949249051953,0.0136592914600201,-0.0364727497643199,-0.0173374543792277,0.00761784669366564,-0.00686234835854841,0.0239635874815116,-0.0766163086875418,0.0628718376079909,-0.0199592075962383,-0.0374363118956966,-0.031733829356606,-0.0536789747204999,0.0283023179826803,0.0735953091800191,-0.0468072014892507,0.028102749114904,0.0559284891076053,0.0303937943048746,-0.0303055777914382,0.024756110135818,0.0731156995115409,-0.016562336206698,0.0482336951644317,0.0417491656363705,-0.0773743706929378,0.0276682764793133,0.0189421737654797,-0.0460986113390413,0.0773306034492937,0.0527068567148907,0.018282513515224,-0.0245420428107223,-0.0158211055099131,-0.0410876092878578,0.0748563408827671,-0.0290534287453878,0.0586476196249237,-0.0694630873154211,0.0225877564133088,0.00982064206610463,0.0614865424397805,0.0661009398038038,0.007565944533593;
//    springCoefb<<-0.96936599719262,-2.7687598101009,1.85929547085318,-2.91329446961438,-0.726320665342978,-3.16187434815971,-3.88541369906735,-4.27705571062538,0.764712907462683,-1.4052480317283,-4.94566432599309,-5.32465393965218,-1.79218638723251,5.57981714840183,2.34284396333158,2.88269330773517,-3.45727014994225,-1.41049764880026,5.77715825770261,-5.81149842161851,0.658161335295611,-1.671708297464,-3.74006800948753,5.7544138552329,-2.11024295353373,-1.74018213733818,0.965474436209481,-1.95490890420149,4.13125865230887,1.44059471837607,-2.17424060355661,-3.12129255073542,4.95502010827518,5.96824006729657,0.248598179650196,-2.0544857570678,-3.47681948671481,2.64636968058284,-0.0483562676931862,3.57107862659955,-5.04206355114546,1.28916460631372,4.52960988694737,-0.551064738377972,0.585796554715555,0.589268644427276,-3.95155663649447,3.4117116047733,5.46195619562702,-4.45758357879186,3.88339838315479,-0.16306767492904,0.153893323744135,-6.1398548321844,-0.691839025547816,4.32683557021041,-1.59685176367091,-6.00954979519001,-3.91125853399108,-3.74877831721371,1.71423012903774,0.197686062452311,-0.586885662097464,0.386065031461244,-0.117259070251117,5.94489771755273,4.61476447439344,2.68910664303407,2.30808219813557,-1.71677699329975,-0.0229999362180522,3.54920384699011,5.85557281301398,-1.77657525512236,-3.28504609723953,0.158184161877859,5.09587858930492,-0.953417533734001,-2.71328943334884,4.27464958493194,0.87218408747414,-5.11307625604571,-2.17160329584877,-5.25710778878172,-4.96974588237845,3.41974288445509,5.35291298142868,-0.283412446049358,3.69337828926508,-4.84153075841407,2.25099443673693,-0.875576787548853,1.63934050988991,-4.61907643121221,5.79367344976406,-4.76110376621288,-4.95736391365948,4.12525271830583,4.2111880768212,3.6339034844761,-3.87470948084559,-2.09499705939685,0.899922125614542,-4.30232186783161,2.41161288548079,3.89806122837501,2.13904749404625,5.67913220915637,5.85612577848375,-2.32027571960238,-0.756830599652277,3.58169676087128,6.03519026046145,2.7377237109806,-0.985226796749711,-0.763932973182515,3.53962288650933,-4.37712329057574,1.15114869209827,0.145080939447743,-0.909782240118577,0.778534680131379,-4.00289343295205,0.333707252789741,-4.48760447129872,4.05295768965864,2.2142998993167,0.309284087483526,5.23113748815393,-5.03676653552853,5.36917454446111,3.9228494853627,-2.48021944713258,3.84753982147958,2.86458560063407,-3.96653239988665,-4.58308977058043,-2.46454052925143,-4.57058539658196,-5.01014919794835,1.49836895114619,0.955769209617436,4.8547327687746,1.25037401160764,-2.58969228525364,-2.41367923382678,-5.79674416742654,-5.33325459874431,-0.507617324402514,1.63758972467173,1.09501154070344,4.86578564133058,-3.86706079519689,3.37530331360306,-1.08369231173135,-2.07148006649561,1.14507579741002,-5.15257761241465,4.52098922683959,0.093028085563951,-3.90615894794318,3.60697856544902,-2.26730763492502,-0.103193195075764,1.17133318692861,-5.68590724014262,2.21345960503758,2.87142861634818,-1.86726256354238,3.9260594143073,4.14446461839983,5.91429158760381,-1.40135658192694,2.71601218132276,0.881480393359781,2.29213633281942,-5.98085225250402,1.36792142593324,3.24206693407512,-0.205284376906537,-3.27767404939503,-1.94610672522145,-1.62268394142763,-0.861549644591923,-4.85398861747006,3.57680895269269,3.35015548891247,2.57427237993996,4.70741654027804,1.58795950990039,-3.61588474034776,-5.48192761351681,-1.08824712465059,0.399992824727221,0.698064397259102,-6.20009914357365,0.997270790436275,-3.37166120355499,2.95451467277452,5.41319342689389,-5.72878699509936,0.815794085322674,5.04429980864603,-0.8469583770263,-2.75137893335457,-0.35740499799419,-4.83800724420688,-2.44904598585859,-5.27266877791263,4.68724488986824,3.62885483723487,-2.26715762730766,-3.54204704120489,-4.27701430419276,3.15447792810042,-2.11285045867495,5.58297984849994,0.221448217012884,-5.82176327873499,4.00721118292631,-4.4737774789384,-3.15446281908274,4.80846876940951,0.721160596411009,3.52871521149615,-0.776652039183064,0.804246652837354,-1.75719920391925,2.13487196311362,-2.5244238802398,-2.62719098287703,2.68927016801426,4.57455541093455,-3.86607638008267,-4.44087340901204,-4.46000872827169,2.05970382192314,-2.99569545321892,-0.625869514130281,3.07022024986219,-4.59163576920235,-3.28019988274708,-5.48012258329714,-1.85049760455556,-1.27402898693984,3.95754055065495,2.31983713676949,-1.97423434429157,-2.10419643818384,2.7812590580345,-4.25020836721693,-0.294788717122245,5.90998144480343,-5.72492479780742,-5.85681332071124,3.15551145044791,-0.218391636990486,1.23061853212612,-4.88487295932301,-4.36670487972854,4.989379857738,-1.22887874220004,4.6057504941374,3.28075006282088,1.18823007771729,-6.11830812072631,5.10392653454918,-3.03525130035956,-2.83081837394523,-1.80512818543277,-6.24821625634905,-1.13926894314758,1.19785713182015,-5.44515363964619,3.29341865229685,6.20701335073199,4.79557211100875,-0.669929410933659,-2.05040619941125,-3.59180953302676,-4.17185555875083,-0.0174293607765088,2.39658695570267,-4.54505931979907,0.540831041416068,2.82295883499143,4.89363733064884,-5.96074580142609,-2.22960783873412,-6.27442083452584,-4.04426547530295,-3.5234131868478,-1.22011437672588,-5.72170018701723,6.04052207597308,6.25130090684308,-5.55682310774354,4.86126340467059,-3.06713559936815,-2.10445628168877,-3.22704998661385,-3.0321666557172,3.03945997516365,4.25399235105798,-2.19413509536339,0.0496934274604999,4.17782049593829,-3.68174818435464,5.66294921652684,-4.15577090932463,-0.990372517381396,-4.79209154807566,2.10998492989886,-4.8769707675304,-3.05396566202305,-3.63236923453674,4.22917326746103,-4.44351353722589,-3.30992983011116,-4.28361977127309,-4.43474916590006,-1.07101010541412,-1.52384775812089,0.62832165737406,-0.509525092431343,-1.76651088214781,0.596437358365473,0.216836999825116,-3.18843268332889,3.81248695899732,4.39556591813635,-0.132297464091069,-5.50286490257155,1.15184068744832,-2.16149031888477,-1.41381479208326,-5.08165109094285,-4.26685502879815,1.1876222235621,-5.70188708026768,-2.13944073812278,-6.08593549967097,-4.21079342249167,6.25372939762775,-4.6797210613497,-0.98157388451472,-3.66182504276066,5.83263740611133,0.858097778259394,-0.688569672871827,-4.73416757101343,2.70653381235934,4.52360542171406,0.0251698767173524,-2.9483297302666,-2.26910487656896,4.54184419456954,3.93129282809887,4.23091732325616,-4.92977369461102,1.46059458709619,2.34329803554083,1.22111404129791,2.24091488452464,-2.78804647701085,5.34280892241314,-5.45608511341029,-1.5865123679537,-5.20723131223668,2.01472231015181,-1.00521424236971,-1.06348684450778,2.21197201633251,1.06717753513863,-1.0929426527317,3.81543615498281,-6.19758155522777,1.52841750450763,3.36488835524247,0.943701423031627,-5.44333736836419,4.91390599008071,-2.63294996460904,5.36345325334986,-1.34410933905361,0.70190550512436,-3.18883682321909,-3.08545035033574,-1.64998687262844,-5.2411047058146,-1.73203884494676,6.09379291446775,3.3853785355779,5.77226040220282,2.0515225989924,-5.68585314728462,4.83188411876429,2.8786226855821,-0.989180315238316,5.90783801237929,-1.38984020426609,4.28879064239198,-1.43883403798017,-5.46105339378525,-0.927217028321069,3.75140850928813,4.63756796119756,-0.841613383548836,-1.00335919205591,1.71927111058836,-6.18109716051721,-0.163511360420103,0.349991894817399,-2.53086192512625,-1.08324330707024,5.28906775576379,4.45422877999811,2.01110506971067,-4.07956780042362,-3.478943298482,3.05318556389606,0.471578560481289,-3.66833558401425,0.155378893622293,-0.039346243167563,4.66637221497815,0.752710946337676,-1.49064733025494,1.26180969470858,6.04671583695103,-1.86599452372733,6.15515469044249,4.05232127349134,2.9783566382925,-5.58908390919443,-3.15808095482973,0.446579941728962,5.33166925200313,2.28349086162143,5.72640595552473,0.767755162591486,2.38557890110422,-0.720290610747049,-5.16543814844279,6.13790217597798,4.47965128218271,-6.159555592679,4.30894575012442,0.207571146041708,-3.95593819310262,-5.45318274835758,-3.02242849006223,2.79882556737867,-2.83833313237183,3.41613560356007,-3.52370588164057,-4.45514612324535,-2.11433865010226,1.26883199395616,3.08984877146323,-2.3508080190029,5.68602267022883,2.96181826190572,-4.58167195136323,2.38119410852133,3.65591955271128,-1.45656770619297,-3.45541115560138,2.70440360471441,-5.45626204457154,-4.01219040592832,-2.81102643854578,3.21250205653269,1.5507041891763,-1.69327938698856,6.02021343710826,-0.252829734492659,-1.56964977966756,1.09297957678341,6.23792661740072,0.757597233081488,1.92298202842583,-3.06768707851318,-2.72676240539152,5.367834096054,-5.93473668080478,0.0327169188195872,-5.37049723304302,-1.76589013090704,-4.98163629307593,4.00253674427188,2.16648705009006,-5.57879882869877,0.681169800325924,3.8680003045785,3.08558047982257,-1.94609585281446,-3.87175260746614,5.91335453007286,-5.52487745395173,-3.04482945203768,-4.38202108170713,-2.0527186924975,-6.11551259550499,3.45186830746917,2.5371871263656,3.23489163700567,-3.08414662702349,-5.31564785915363,-1.95531399206259,-3.12940521547444,1.72513457392786,6.25085324221491,0.0860929060123809,5.28155737438801,5.33550213241724,0.434541425207597,-0.968910912644073,6.2481901052259,4.95183650015223,0.332637994280002,3.9675416436461,0.835138344390612,1.03702436558123,-1.63447376187964,-1.58004655688256,-2.1605803545962,2.70261558530589,0.831386041502963,-2.53041103037501,3.46092333720584,4.06974178946528,-0.629246912082139,-4.87498056114334,4.23741439396029,-3.46056380461297,3.94539176522227,1.18912082511429,-0.261525231636459,4.91292911192031,5.51699203890337,2.8922547528891,0.354878479996494,5.48466007526661,-3.30483754695019,-0.646749351467167,4.53697700183218,3.41288907825741,4.66752493588876,4.50198190120641,2.08154037255796,-1.28302227568291,2.18633833900084,-3.36650648305143,6.03718729575,-5.73132062873048,1.33663216006601,-2.40657826469787,3.25448016242709,-4.1151670042827,1.34619590492712,0.43221829378125,6.23775999103426,-5.56623621300669,1.84042293263791,4.19198917914288,-2.74361481761966,-0.497370502139822,-0.902075195742831,3.27804515074388,-1.86762659607119,-1.66826836269113,-0.112885302218691,4.77043708392531,-2.46679349327619,2.86546235083112,-2.15949746754186,-4.21300169729568,-0.00483377091147422,-3.77515773750477,-5.99420499608927,-4.20647860420519,1.22500519266399,2.47531854291157,-1.28979988725661,0.979007282562317,3.02718312003276,-6.2363529271906,4.85561421786444,-0.00152192339182384,-4.0683347314733,-0.0813750830601098,-5.85248883546225,-4.11375994629072,0.635573903933198,2.27111929717567,-6.20495596714784,4.17514429216521,-4.50943641081583,-0.823845962890671,1.17000423705742,-0.093877806887014,3.7910708744182,-5.22606626516127,-1.60662592296171,-4.95890782470967,3.92258128566984,2.51706180949643,-2.88872432200535,-2.3654376910933,5.02508927199167,-2.59974411809462,-0.288731089446815,-0.0330907411960161,6.15875963066862,4.70465422329657,-5.33726866448537,2.90275754484971,4.75148649610597,5.80153075337907,-3.38194958439378,-5.60003343536733,-0.56302953553271,-2.95125321985603,-3.43060818165805,-6.21064083159951,5.60305127731964,-3.35237894880589,4.2476886605657,-5.18957033349619,2.10696028830344,-0.865492308228559,0.999737059616794,-0.385154043130037,0.191626626610166,5.67629633665509,0.939123332160296,-2.16897728771999,1.91017294029985,4.33358421600662,1.74877022118671,0.652077012291515,-4.54934510793967,-4.82314607411178,-5.66419893475617,-4.67377068312272,6.16469335503646,-4.71828239338987,4.51217206172699,4.63299464529076,-5.19993684586247,-5.15296272851846,5.31614640992343,0.520218824456492,-1.82103074837449,-4.39764697173462,0.592763192856979,-2.50116467105486,-1.46684072054052,-1.442733352429,-1.40754980455105,-5.64306563223708,3.97495953934244,5.87537245506574,0.254965524632883,-2.11659903404739,5.26848359172083,-5.08909634320682,1.99760887823262,0.895471326169011,5.5276730727998,-2.53680610643235,-4.73563686739115,-5.30485724099154,-1.07676697469245,-4.11665059629565,-3.69544271826258,-1.19525882550766,-2.55174778968552,-5.46645586238726,-2.84544938606857,-1.46849943554799,-4.33623338505405,-3.81248817614514,5.3349045889085,0.125921066571457,-1.92694994787977,-0.355517424086193,3.9079415955166,2.88939453157972,4.48493442348481,-3.78279341488612,3.52951409934264,2.17670876282725,-4.19060615982038,-2.49870557602448,-6.22307547707181,-5.20530777395121,-1.3046167192313,2.05771860116081,1.9733487522178,-2.06012885228317,5.80409769472846,3.52089709067832,-1.0818008932747,-1.55585448581566,5.68743169438268,1.50594158846271,3.53207188867668,-3.14750130115451,2.32267093192713,-5.59656269739189,1.6671844632975,4.26962274687308,-3.12586567353704,0.718903852205999,-1.88764139240714,1.23036957858319,-5.91979877773187,-4.26288499689054,-2.16342109568876,4.84832085160462,-1.76249321177666,-4.91709219634613,0.741844408580199,0.330085828402965,-1.13261257237061,0.80195413150839,1.40796325445175,3.8459559083981,-3.42351247318247,-2.90187319333045,-4.49735814973674,-3.90259998430568,-5.6641613085038,0.704026156988551,0.824730735730331,-6.25991481997279,-4.07321745454873,-1.92638258144466,-3.12423092112731,4.53263867737839,-1.23976007883656,4.82613874217019,2.5190762183998,1.9175594476264,-0.738142611475482,-5.65175037985901,-3.13525617964208,-0.374756183355675,-3.63145017674955,0.984507930520832,-1.80962053760273,0.889241817325465,2.35060093417471,5.21540907097747,-5.06385756012324,-5.06519684404757,-0.265822003365814,2.62729090018018,5.06394426435052,2.59385072345171,6.00860290684973,-5.71659909123789,4.97443593914604,-5.93874360750574,1.27061226575066,-0.484018530975306,-5.91547322747853,3.48058001120193,3.8727840934317,-2.75651894860584,1.73003348272865,2.26026862915325,-2.10288050243167,-5.24006148139457,2.20164648322679,-1.19121488697991,-1.56711121307677,0.896505386870348,-5.13517681657608,-2.64171509080058,-2.52189299689462,-4.70350557514983,-4.63841754793216,1.89292590300054,6.00760090160882,-4.17304122457604,2.34127251431192,6.20475875409068,-4.51355148227108,-4.72278108292085,1.71028652328188,1.39932155173067,1.92099004313041,3.94325032236668,0.307569424637958,6.1472880300036,-0.942256428648643,-0.58885120368612;
//
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
