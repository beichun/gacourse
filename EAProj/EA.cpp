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

//    InitialHeight = 0.5;
//
//    int num_frames = 1000;
//    int skip_frames = 32;
//    ArrayX3dRowMajor best_position_history = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);
//    ArrayXdRowMajor springCoefa(num_cubes);
//    ArrayXdRowMajor springCoefb(num_cubes);
//
//
//    springCoefa<<0.00777953353141413,-0.0561735179909382,-0.0155445242373014,0.0584407754142027,-0.0552563985321933,0.0716803874036672,0.0541639547674749,-0.0685322850889211,0.039785556023151,0.0678447089287661,-0.0428364026373422,-0.0427289237653506,-0.0691219451972851,0.0287276429444215,-0.0423353988501874,0.0143773018263175,0.0796877884863353,-0.0429808368361466,-0.0756107677498883,-0.00726712521504012,0.0728016649153091,0.069108276213104,-0.0424121904943195,-0.0761236993950436,0.0089584295120828,0.0651174053480464,0.0785027006261529,0.0716475664785353,-0.0689282268234194,0.0780684234937972,0.00612183716433208,-0.0745955240701304,0.078176278465666,-0.0677802150639613,-0.0401988794469269,-0.0339367062383968,0.0731432416258116,0.00361623147444086,0.0489290588949477,0.00386423748166498,0.0647722218114753,-0.0769645383939913,0.0119348110102171,-0.0271124670408259,0.00254033636420049,0.0286765730227089,-0.071478875498976,-0.0748556133335715,0.0424993251089469,0.0657832708888609,-0.020621644975907,-0.0669325552819914,0.00235027399023542,0.00553928684701178,-0.00891205419269952,0.0108821873790036,-0.0647280935126953,0.0315951024515532,-0.0581696980344922,-0.0357682553472781,0.0193914493636188,0.056228277076142,-0.0257687331483554,-0.0478739432095894,-0.0347046845195371,0.0103192393902313,0.03819513996979,0.0204746253511285,0.015306479677235,-0.0399050301126694,-0.0157627023457376,-0.0702706312342876,-0.0791860450800443,-0.0483315469922179,0.00925473388715401,-0.0512633004650768,0.0716563966086397,0.0743523391309904,0.0604176021276124,0.0171817057287235,-0.0774569921183665,-0.0751686205692443,0.00352792184964191,0.0148762261564267,0.0662883771147059,-0.0548966885241199,0.0189873311384522,-0.0234001271908172,-0.0535987535741174,0.0366775375402894,0.023848188940365,0.0583746626127393,-0.063014188028413,-0.0549912978592288,-0.0151818914409643,-0.0669404713748677,-0.0559938552863867,-0.0543987479686731,-0.0770114544765146,0.0112688746728324,0.0112178256042385,0.0770593226687328,-0.0553158394877407,0.0796665241567728,0.0607076157539653,0.00261047545942034,0.0228758435430358,-0.0668720465837382,-0.0531870367998197,-0.012141789091817,0.025529787962106,-0.0639967232123002,0.0345350880445593,-0.0116624820659228,-0.0763355804217679,-0.0465349987738463,0.0219988228110591,-0.0431719595767427,0.0195969851033748,0.0725901858474082,-0.0683039283697977,0.0377214984398901,-0.0589445269754829,-0.00631534707095257,0.020369018226091,-0.0618589278226061,0.0246391200016435,0.0712791423002576,-0.0642309976295712,-0.0706126861975587,0.0248542460177831,0.00416637007341086,-0.0787714038783551,0.013958003266695,-0.00892107483414983,-0.0163063040824171,0.00286160361155011,-0.0792580787834051,-0.00317835059164946,0.0296745668117304,-0.0113998678752221,-0.0622090739673977,-0.0253802921368649,-0.0500804475555571,0.0725961185957287,0.0410589492884739,0.0595514986568836,-0.00600206341873951,0.0382223309381969,-0.0687116825528032,-0.0264085175592492,-0.0120598452221881,0.0138859893818786,-0.0492317944621815,0.026440697883461,0.0415757810145504,0.0719375415999248,0.0390192756611012,-0.0553230681900508,-0.0475555062887983,-0.0338494866126447,-0.0215040784801841,-0.0247658312808563,0.0320717087909913,-0.00723788123914873,0.0268371425321499,-0.0515486852692201,0.076730257569221,-0.0240098853893624,0.0248746543914358,0.0647526769269038,-0.0368893127687691,0.0447366539969745,-0.0406276152099611,-0.00696976032432623,0.0373327725181974,-0.079568665995993,-0.0274182617419484,-0.0486692909750479,0.0386536650167097,-0.0161299442202458,0.00492219146570292,-0.0534061802799842,0.0777560451616328,0.0356903970035214,0.0530345176779826,0.0393318261016774,0.0733445596663955,0.012053793264578,0.0640087579861324,-0.0542109466224028,0.0582043067264391,-0.0374953205685575,0.00102322209674083,0.017348373428615,-0.0758310777674574,-0.0592178265653634,0.0708912982004189,0.00944629777662749,-0.0125922104402409,0.0589083389094604,0.0431170251048715,0.0633601643812657,-0.0142510990492306,0.0488598865125607,0.00317087438105181,0.0356862956079125,-0.024588116083568,-0.0784781321131057,0.0169218028043033,0.0414867269006171,0.0660672466392011,0.0692274216119426,0.0245410715530352,-0.0264118661481942,-0.0178249347665463,-0.0291007186235398,0.00226939063624917,0.00944236162926381,0.0135320608753395,-0.0052390217432934,0.0641896746233989,0.0364613906836423,-0.0384974801347114,0.0236065238451616,-0.0261902358877427,-0.0343285579021687,0.0443886972797982,-0.0352989377618297,0.0551177399489646,-0.0482035132349485,-0.0563905989268751,0.0182347649793302,-0.0648433488536828,0.0701796694799232,0.0401802124642675,0.038506033252229,0.0118221290327991,0.016733663872226,-0.0289208224177923,0.076995347662361,-0.0480592197776116,0.00591663997895859,0.0372345933864054,0.0273366404079537,0.0595047738307644,-0.0605903413801409,0.0782359218589198,-0.0182258355329865,0.0287146088242133,0.0117679826597534,0.0565351427237201,0.0129042833731064,-0.0317706266566043,-0.0505918606978803,-0.0448234628722181,0.0445424920900457,0.0196466096209579,-0.0682179421131583,-0.00616177533108823,-0.0145635116540657,-0.00641015786976095,0.0162923029513528,0.0536852117505321,0.0101322789164876,-0.0476428261248594,0.0485162695350667,0.0276096264028967,-0.0685884788579254,0.0740081769945138,-0.00397378016448291,-0.0184785581838705,-0.0226823507727507,-0.0317489810808324,0.055588326070266,-0.053345525252328,-0.0471768754567843,-0.0403329052591384,0.00144908126511103,0.0140074941953679,-0.0143953321568646,0.0156288189141214,0.0658142154597744,-0.0665116895672454,0.0434549950638111,0.0236837628221529,0.00141042318260782,0.0523391423990667,-0.0501886919001996,-0.0092860018132655,0.0171052226131341,-0.053749604957993,0.0567638844143431,-0.0229357587280384,-0.016233506620039,-0.0224541884392752,0.0399583457456322,0.0758489860574943,-0.0758411896954482,-0.0659438244187943,0.0437114810774622,0.0156336661128484,0.0465786204331455,0.0292478792878091,0.063884685032016,0.0221669464289057,0.055902354035481,-0.0632921904247683,0.0535590057900875,-0.0226485647739137,0.0307153037705996,-0.0325612909870973,0.0729802541402077,0.0165295191558681,-0.0190729805543427,0.0364352491295129,-0.0397867180964848,0.0504061727875827,0.00877439145407378,-0.00997540992217856,-0.0269485591850004,-0.0541203859327922,0.0162749851198285,-0.0501846748451631,0.00294385533916943,-0.0799585215002105,0.00736113671556168,0.057590904895957,0.0758904645572838,0.0115199470946193,0.0716470804771628,0.0396019455602402,0.0428881060686591,0.0382257008358025,-0.0111501752264566,-0.0689617018350222,0.00949843863729949,-0.0352478212654813,0.0456269476269755,-0.0377733116400304,0.0221036139606049,0.0584614115108091,0.00966539737287229,0.0150838680263068,-0.00500906940782865,0.0705924168185296,0.0378133369425447,0.0352042125701924,0.0529298594467947,0.0602935086098935,-0.054771197426492,-0.0540186997382057,0.00569924204046968,0.0415037876933365,-0.0242033745833688,0.00911697794176497,0.041545266193126,0.0631577622066987,-0.0132921172367838,0.0546539668307918,-0.0053222907731879,-0.021645036759621,-0.00296232376385588,0.0218313223597739,-0.0634193359983244,0.0658875010096875,0.0328696205247518,-0.00302668873361624,-0.0493603202557938,0.0606157282649613,0.0391999996263534,0.0117082900339522,0.0688987688668532,0.0142428195170326,-0.038420360217998,-0.0697590224397178,-0.070658341045798,0.0357710461857594,0.0207623119376424,-0.0189598383470251,-0.0774997510004322,0.0595228049808754,-0.0439679476450979,0.0133010580080101,-0.0725052431190877,0.0776746912289759,-0.0245252018163564,-0.0479922556728089,-0.0743984237287186,0.0334503016404111,0.0654005583307708,0.05469849090578,-0.0034726472587663,-0.0610183037729088,-0.0353583166540406,-0.0103470320675275,-0.0767296573504478,-0.010029379171333,0.00573627073584882,-0.0733648889480926,0.0547580005856035,0.0118126780054898,-0.0142485793373774,0.0121325891893975,0.0343716472077551,-0.0481453252435407,-0.032899518214585,0.00601199992402206,0.019507547067249,-0.0250667284918328,-0.0246442332606037,-0.0178975908914104,0.0667140048680427,0.044402779569069,0.0743016018691946,0.0282553133089996,0.0208346737767096,-0.0644519603925068,-0.0253194889916663,0.00459802504843008,0.0488634010200496,-0.0399189307354013,0.0703010992660658,-0.0344743060853678,-0.0209372344338043,-0.0450572174624806,0.0351786618471046,-0.0176668917842521,0.0249134033661864,-0.0390850674915524,-0.0110317807323447,-0.000328596122715901,-0.0227281319828369,0.0547196399302779,-0.0681960070078243,0.0718686459920689,-0.0419906017752321,-0.0469639398748819,0.0261796825128513,0.0271960716448706,-0.0142592309109211,0.0415308806866086,-0.0611170060099647,-0.00203171704059081,-0.0670858687847321,0.0702168231280454,0.0176905537851576,-0.0180279271016027,0.00433375295453415,0.0266215058912623,0.0396212709460993,0.0520345107335758,0.0759893219154278,0.0248522633057331,0.0459851571014082,0.0614978216773419,0.0651923133363911,0.0619205025499316,0.050079764337316,-0.00821809370453381,-0.0619067872231392,0.0628375478008937,0.0111741499282299,-0.0150075187976507,0.000862329382560935,-0.047732187047476,-0.0231388728800876,0.0267570544971443,0.0635126542968269,-0.0302100624657283,-0.0642739925460303,0.0173504005825847,0.0202250982646125,-0.0560897216834546,-0.0420978483381206,-0.00364643808624075,-0.0536038489516842,-0.0788733305031775,-0.0652689655801603,-0.0355699366869265,-0.0609683016599008,0.0156840563265998,0.0482188877315348,0.0580806850539896,-0.0161076248139644,0.0319181120171762,0.06659123934181,0.0790812389611681,-0.0125263160898007,0.0139408734133192,0.00685922047812542,-0.0470513254995697,0.0663393522735403,0.0433420395579988,-0.0572496259106554,-0.0690748070129542,-0.0788533683860923,-0.000388498790743066,0.0363659099286264,0.047819905336552,0.00265231076751478,-0.0709973016339342,-0.0425048531417292,-0.0731216936526456,0.027961225932446,0.0432431607522272,-0.00822675329084823,-0.0212925228948204,-0.000608118586525377,0.0464557718329391,0.0792163840677665,0.0374512196879141,0.0795779152212562,0.0608981437147121,-0.0389821964311331,0.0282250730452244,-0.0361903691087805,-0.000839072801563461,0.0183835252921482,-0.00423392888355718,-0.0603755542544534,0.0266308570404681,-0.0426318783884085,0.0635474534121328,-0.0676134914102096,-0.0776253193605809,0.00545429032549928,0.0174321090883725,-0.0103329117271737,0.0429952481216729,0.0631683037351669,0.0799537194241347,0.0184131503563436,0.0681520107293017,-0.0128975988612033,0.0282158963141106,-0.0293193711849485,0.00968279581159818,0.00105755548973454,-0.0580664847689059,-0.0123957688791658,0.0592271626165573,-0.00596409451494184,-0.0789510233695391,0.0693176935303129,-0.025462382261391,0.0108584170466561,0.0774787478882255,0.0440657045524873,0.0397039381972067,-0.0658402243190632,0.0223101222246467,0.0769529631719705,-0.0565366187582429,0.0714922543575486,0.0212097762996376,-0.0100830512168273,-0.0652735031327575,-0.0383569707527556,-0.0127848328337003,0.00383254852324377,0.0734338938414277,-0.011359960404858,-0.0554434634630771,0.010685449764995,0.0709899546164041,-0.0337639930629004,0.015242421447878,0.0651131162723122,0.0466229545169617,0.0471601470406913,0.0256920496866536,0.0505474099379719,0.0707165940807744,0.0165011487605521,0.0378347295326343,0.0715724837740755,-0.0736992263019547,-0.0723826581949287,-0.00107631404002957,0.00291990939663719,0.0587947419932088,0.0763286140544008,-0.0623483247972784,0.0148078055799425,0.0391318833637665,-0.0605681576489323,-0.0236640476918146,-0.0372392857620676,0.0365560127592441,-0.0523759519925229,0.0398466269484938,0.0671201171666012,0.0661770662600999,0.0751023334241949,-0.0546794636429657,-0.0303098696983931,-0.0266636703660077,-0.0783488497130335,0.0382834726377779,0.00651328753983289,0.0156465914545798,-0.0271963482849283,0.0201049855817598,-0.00364833410999195,0.0593428044581221,-0.0558802953995207,-0.0790248783673276,-0.0337999557488598,0.0178858986595358,0.0368069238759609,0.0383390541925743,-0.054094766217328,0.0144425783746143,0.0224534506827842,-0.0738687648595631,-0.011204484743627,0.0495103698826909,-0.0194886500106606,0.0314345144254316,-0.0511480313405153,-0.0199417216051099,-0.0265158752848003,-0.0760288450103387,-0.0731350709465961,0.000755630024129346,0.00556260790934907,0.00599501800071216,0.00818572316699928,-0.0623494353808227,0.0222753911196605,-0.0369942145594369,-0.0689062280156213,-0.0397749040833557,0.0480031458486556,0.0659616766804651,0.0313523989689315,-0.0459677964290454,-0.0179571804115349,-0.0587858902564207,-0.0222342959522429,0.0682598848718376,0.011602690635064,-0.0440740962904292,-0.00173710737458295,-0.0459438587566576,-0.0379428611499923,-0.00697752166864347,0.027279850368984,-0.058645604764412,0.0206485956630896,-0.0343747230406733,-0.0428155311768947,-0.0270761027685721,-0.0473063982125867,-0.00122839491871576,0.0711322095203829,-0.0349175370274659,0.0489970665816967,0.0626768513688244,-0.0208775255553785,-0.0203903893103778,-0.019265219317407,-0.0711909227032172,-0.078025998546754,-0.00110230442187856,-0.045152653886542,0.00903759364459551,0.0261448852467094,-0.0202557651420384,0.0152434243705326,0.0655057451806524,0.0140385768255399,0.0417620502979318,0.00217261333119712,0.0292024722458806,0.0470615212465923,0.0568330892812615,0.00610163723645652,0.0460887782862823,0.0378716036481185,-0.0706798221127502,0.0342766892697088,0.0701652163134808,-0.0553692505393966,0.0690015242197558,0.000262560521419878,0.0406523964370845,-0.0383345756485754,-0.0292788603665674,-0.00361785132606414,-0.00546493713067143,0.0499382525542463,-0.0137306160916251,0.00468369025955148,0.0473534002003043,0.0245727538059339,-0.0119739135783044,-0.0059982251031316,0.025160766665433,0.0232808883224059,0.00946580436521479,-0.00363824825903319,0.0689619390615085,0.00336972501285826,0.0034886405260715,0.0329058823515223,0.0677110256946599,0.0526726682356897,0.0698723039728926,0.00656969916381394,-0.0367889088936098,0.0474873211828467,0.0274348883831105,0.0709339829492075,0.0780262356987345,-0.0790063206846855,-0.0625715202757025,-0.0603026404140064,0.0581505505661345,-0.0260281518595424,-0.0559546651579229,0.0297056130154052,0.0701248273207456,0.00912809078075368,0.0239714673645662,0.00825180795427961,-0.00689145770244834,0.0653916191241665,0.0354163822324091,-0.0617306911115212,0.00867250737206662,-0.0351178134768819,0.0146310606294456,-0.0023655536409307,0.0482519115359764,-0.0618802989189887,-0.0494596713639142,-0.000418301634685276,0.0707923693912068,-0.0595873674655274,-0.0738486024708713,-0.0459965395769088,0.0678999537918251,0.0335862859122391,0.0195898195673648,0.0659261894160538,0.0345799653020594,-0.0376340769779096,-0.0743764509979526,-0.025785101477888,0.0163377712370538,-0.0503311161558754,0.0716135602405358,0.0564184381637606,0.038796974699384,0.0155850275305961,-0.0652855936369326,0.0415878329855125,0.000976646580256819,0.0501307885954766,0.0157090360933428,-0.0703508460476766,-0.0649870249559111,0.0648058864403543,0.00728360031139273,0.0632648865800653,-0.0770744125531402,0.0378239290219843,-0.01715341505462,0.0737179568380667,0.0582365615564569,-0.0110020175254913,-0.0522785827388422,0.0461365152737761,-0.057415731687758,-0.0273411394410493,0.0320627046898299,0.0571642336143014,0.015024783655547,0.0376862536918774,-0.0486208678635866,-0.048637445181905,0.0673551375360019,-0.0570073076975566,0.0378251533013886,0.0261521121608802,0.0385777198330396,0.0525395597389618,0.0580576290833101,-0.0404456336612094,0.0226703482599325,-0.051767545180287,-0.030796479708886,0.0376833233785272,-0.0669616588144385,0.0564871206770125,0.0209482098840867,-0.0640360712930728,0.014311049624491;
//    springCoefb<<-6.12712579493439,-3.7355741534378,4.11994380623928,5.23167732247124,-5.51412339118599,3.86597615264088,-3.96189739565131,-3.30505102765612,3.70182562715742,-2.65219548255993,3.84630314077596,2.68557631748042,-1.18012699801737,-5.80180280932392,3.66664053597396,-5.64817556404147,-5.00387224374645,-0.660654873268676,0.308307460945029,1.76312520407976,1.22614366475275,4.11696300537213,-5.49163718204964,0.657611529728839,4.5531938463855,-2.68592203427518,-3.42391846220109,0.177590156087718,0.116649811694589,-2.48441985937549,-3.4707366617388,-3.48247265850974,5.91600368563545,5.19301675192188,4.30152083827309,0.454628770575732,-4.74134797074105,1.6727155495837,6.13984848312005,-4.10091164025752,-6.24966406538637,4.87172333898386,3.45023133129733,3.92276722255463,4.67879137692575,0.611248919834527,2.43551315706077,4.96759306972374,-2.52078844694973,-0.569218866794794,-3.65079768349557,0.429684168609669,4.26789642749645,-1.28803581874249,-0.583509920168844,4.1326160296037,-4.62633863146579,-4.07797314952361,1.62024629862761,-5.12175710060228,0.525009126220878,-5.92470568576858,-1.11354079641959,6.02872489411741,4.57415720463998,2.79072538425649,5.50237244339253,5.8214437715801,-1.93945302280701,-4.34074544565605,-5.85976288535207,-3.87012473662294,4.87327549504021,5.69912945632733,5.14326175024045,0.264736921934392,0.604368077060155,-6.2726839645842,-5.92069099389859,2.31323689523245,1.86803382564432,-4.05980850565181,2.6880951288258,0.16466377704376,-5.50653218066184,-2.71075983151988,3.6308153931988,5.91630357391487,2.45779317082311,-5.69630730599741,0.341134891520327,1.29671652871059,-5.62845075317405,-3.65321634347091,2.34183277681695,4.52922805363167,3.81599551633032,1.81645552482255,-4.74068190415061,-5.69277287256986,4.20522360191584,3.894226346689,-3.66988981156256,-0.443997271084643,5.25040423597164,0.696604340671824,-1.52231006154659,-6.03957542169969,-5.28252666878122,2.8706780047975,5.24975704102411,1.52129911495053,3.49345478945735,-0.0141936903258181,3.761477785964,-1.81902457862561,3.50884039301137,-4.34246224436198,3.74263501520531,-3.82928576405626,0.092009999927775,6.1534856355302,5.42051859699224,-6.06375904295191,1.81271797573654,5.49956565370219,-0.367542938191885,2.3536881839675,-0.941487425292183,-1.39207308535179,5.66917234717609,-3.99803796111415,-1.50969104158599,1.44853371718695,-0.467151817260255,0.352807011437311,0.0929523824804415,-0.122519284701458,0.596416795589292,1.0936109195509,-3.53502648575563,-3.05081979718592,4.10347311614506,-2.41604220928927,-2.11205266509117,0.00173106800863823,5.91940521766504,-2.79173127500776,-4.69123594098246,-4.34523439215278,-2.88464357645934,-2.07199641220063,1.97237692405052,4.41740475092333,-4.34844123141445,4.90914275815201,3.04973084097901,-0.663227999835604,0.972603294408188,1.96307840038111,-4.04446855173568,-5.49350929583963,5.9347011201927,-1.50237016661506,0.490235862488962,4.23217103817762,-6.18735898814643,5.21754490793107,3.40341027627024,0.750067284694982,-3.86745132994918,3.81091304642366,4.58685405769278,-6.04716341965579,-4.88831436286561,-3.80838380739838,0.237752848352844,-5.25209435105224,-0.316929882406146,1.82970211322206,-3.31414353735336,3.08161174113452,6.04089090102143,4.94141858669717,1.21583129205785,-4.59073553624469,3.56737614484918,1.22622908746382,1.02922166391971,-1.74320576659431,-3.09387771800675,3.26793831218403,-0.953529862433935,-3.44236180366571,2.78486752681898,-5.1962446864619,-3.56054691567482,-3.71063486479011,-6.10435110420927,-5.75685930416366,-5.06800945615177,-2.65229444018885,-1.16083343198588,-5.77741865126473,-1.78192173492456,-2.48132225327916,-2.93501201479126,2.69946755535004,-0.262133442743659,0.152822251216009,6.18490474844877,4.71197471205714,3.81849236859271,-4.33868397175302,-5.81036563150348,-6.18106776171598,-6.21524927509465,1.78367268398957,0.151973303900076,6.12133808734346,-3.10720443920563,5.78051623395437,-5.26236790960531,-3.10213359072,-4.96643009668339,3.80568482306535,-2.0151930771819,-2.24379180650654,-6.18813524757644,-1.8363589755395,-1.7174659106702,-4.97295950372821,1.79453178427165,3.40488586319559,-1.54956657728336,0.379031994395138,6.20146504675491,5.52177564266097,-4.54357901057424,4.06518888419119,-0.160975146201538,-1.35330270817323,-0.909780332473433,-4.81216707306591,5.70650516731864,-0.436960763976913,-4.71004963478188,5.77444109222399,-4.93647328583902,-4.18991148440445,5.61259397371577,-1.76049252504465,-4.69258045630175,-5.93295913588954,2.21641360524597,-3.45679725290347,-1.02949759367157,-1.80678831654492,-5.96638396953345,-0.41666223761003,5.29090114318703,2.99207062901058,4.92158149394319,5.95171092484855,2.86419995787976,-6.05783259221263,0.722870964129169,1.20591734216939,-3.88287180872673,-4.8017856740455,4.40415844340409,-2.24759008013115,-3.83298598120948,4.59919683191123,-6.06541281914803,5.19722605870657,3.54089906600044,3.04734806442621,-2.77836845484295,-2.24025011115377,2.43111582434432,-0.138080244573998,-1.0552309776462,4.14538181535526,-0.32174123263551,-6.03061128598762,-0.921235004250722,-0.887925770391521,3.61542981379436,-2.7211952999916,4.01663501430464,-2.039090311744,-2.02454549475246,0.910794011800428,4.22047212930521,3.58918397547674,5.31595469050772,-5.10999048203061,-1.58163904851553,-5.08112606135626,0.29696560355647,-5.68976089399312,-0.260247717964083,5.34571211255164,5.19262113791811,-0.0424753371121114,4.25975297125821,2.45033499806688,5.99185520032731,-4.80180068943642,-6.07310031308689,5.43598855165841,1.34330426598959,-0.845146090733085,3.29818516116199,-5.2616221724976,-0.592572176720702,4.73220033197337,0.133637257110882,-3.26032756877801,-0.344245143080327,-2.13291292858448,0.983767325329655,3.91439456216721,5.06106628321595,-1.07894575121681,1.22039333179228,0.142390334709326,0.0942489667525828,5.92193948327675,1.34444947335307,-5.89197063554262,2.49626218382508,-5.19898345046269,5.73692668286069,5.42479972135007,0.462505159192856,3.71349444826723,4.77115869847055,4.0465991397393,5.19487895883081,1.80203440633007,3.19940248554603,0.254998018968721,-5.32629689025469,0.214402440856358,4.13193001218058,0.364316133024608,2.59135259776763,-4.87298690226967,3.38717377009827,-4.0360777453127,5.23986497756151,-1.91224411042375,6.16150201685451,-1.94483354763819,3.29199533835944,1.09871014864679,4.48074198707114,1.33728271734889,0.737464426071875,-0.457993745427468,-2.50572633043059,0.969643015355514,0.626208004109845,-3.05198485342158,0.111257530853917,-4.61525078931663,-5.62167561100602,-4.5799781555808,5.71453355042267,3.5923629757176,3.56136772973712,1.66193261088676,-1.73415636616631,5.78647463482397,6.04317326384742,4.78795611280108,-1.1077108865686,5.38637283147916,-2.40055727434221,-4.01690567755151,2.89120049806693,1.27418090805401,5.19388994430222,1.41503356081822,5.47631003159232,-1.20581554062778,1.08463257597867,-3.70477495926732,-1.19160856993186,3.76296266074436,4.57964283135591,-0.322566464710331,-0.473157382607943,2.04017536391707,3.54180552513735,-4.48372998549512,-5.81156507801677,-2.30903368003854,3.96920646779938,5.91994224539764,2.42990143851185,-1.53933392191911,5.346059837151,-0.589540776443867,1.35424714466284,0.452971219123134,5.76670258340672,2.50903136660378,-5.48381579014704,1.29617689920181,-0.0590231346911729,4.28137952039722,3.18927104655618,4.60565646957841,5.98212762388251,2.19870998537002,-0.208647994170317,-3.4627537929413,3.6128867493634,4.88292863589782,-2.09562476242422,3.0735778574586,-1.72282302881251,3.71440306081951,-1.16943198447601,-4.46420270952683,5.51385827532439,-0.697811862492777,-0.490051189565372,3.1998795372721,1.44234295811684,-4.34333495105352,-4.62263959049869,1.5857987916365,0.202068836208283,4.36981658250325,-0.229104111921755,0.762901036638195,-2.67763516525417,-2.94081938985253,5.20636116272589,5.10756301528858,3.52185103028363,5.99019512682983,-0.762904465718116,5.11636981123801,3.26839764504661,3.97311657273178,3.73019565653917,-2.87502082918723,6.17772198620079,3.62753273382018,2.06079522901223,5.50991233946513,-0.061535550158819,-4.16216512518616,-1.26258033653326,5.84589459188923,2.23534821214847,5.6595799465477,2.74892671989072,-0.887147234661015,5.56981059866114,0.492113757731108,-5.58453364302451,1.26888550338329,-6.09346435477594,-1.13350255905021,2.29997360030875,-0.56964023400801,5.82496052019952,4.42997734956821,-1.38802444793125,-4.80460490439098,3.46606693953218,-1.2092411212487,-3.59361432332786,-1.58231869098551,5.04407941711521,5.75201481591451,0.321516961446006,-2.53149188689909,5.732146598116,0.467849629172488,0.600301919732989,0.672945766096758,5.1793297666711,4.42927277312289,5.25200295827171,-2.32688792195744,-1.98691717088606,6.03679423105411,-4.59161222207334,2.78545838771805,-2.57948598582107,-3.89296066509785,0.509929501498465,3.29285903162227,2.85977197995589,5.28048685191707,-1.81644568230988,4.20600282814367,0.898170995147697,-3.92051740672953,0.114224474961249,-2.62289415848262,-3.0196377690638,1.81704376417246,-2.72167097072302,1.46130825704014,-4.86186637267948,2.11179193213206,3.4751995156291,-2.61312651210301,0.184401209968158,3.03204686613518,1.44930679162514,-3.39777085005336,5.04798358845198,-0.447906053539125,-2.16427789168354,1.8797114357919,-0.526180865066207,3.6809633346698,2.96723006392492,2.46823129010235,-1.14426102227346,-4.76035925210466,-5.54496285071301,1.46517942454699,1.31032934401598,1.92251425957599,1.03232228947891,-1.7370542722448,-0.949739389088397,2.68039741378677,1.41956704439222,5.95272719061466,4.97752286008889,2.64524169555052,-3.33881319264849,-1.63244591753686,1.03824666175183,-0.953376870463744,3.89354316980899,-3.04121050290825,1.25932276327345,-1.69207784557662,-0.0354939644600374,-3.34834772596188,-4.01121479172912,5.11896258772025,-3.57011297194646,2.38809476320681,2.05730915965558,4.91553172506355,6.06927646180301,-3.78471633289433,5.57163081110957,-3.07459889818925,-0.424778068103719,3.64611759797744,5.81558755290426,3.42695630516904,-5.18038976783539,-3.56647372921778,1.10972994755138,2.2197643195725,-5.65071501423477,-5.67363785174549,0.479074382896583,-3.1886276317889,-5.3721051613873,-2.54257481949526,3.58279245138444,5.13184059684164,-1.50816505420781,1.4789114622208,5.3239654904155,-3.58783900957869,2.84979336323773,-1.27427795803778,5.46899905805128,5.71522570961954,-3.34862935524248,5.45038695395934,4.3537907774929,1.67164258903665,3.31301461519366,-4.82491559053274,-4.43364749497789,-1.23189486961238,-0.20367884774787,3.68520374894401,3.17300838340938,-2.08712051548358,-0.416781693646729,4.17077741740266,-2.44635941292233,-5.06445003608384,-4.31064616999981,-2.32913303888975,4.78916315478491,3.38870103361652,-1.62497961095625,3.76599140311081,1.07861231708204,0.674892064404347,-4.60866366137703,1.80952536015407,-3.21357049580832,1.96596460311515,2.97003641473984,-2.78558255435089,-3.7159379134811,4.97343611657886,2.13469272639352,3.49689075753936,0.77422500546478,1.63559154245543,4.1092159355584,-1.38155098287285,4.30495957345316,0.776266858710856,-5.80642563622938,4.48753907886166,-3.82714595966889,-5.14475987195552,-2.75734366573351,0.932197790312874,5.93972182569908,4.92948744250937,3.60968222576761,4.88257556940276,0.335845957062731,-1.37557068504581,-3.22887052376489,-4.60646404119623,-4.2083662366299,4.61951773654755,5.04027591794731,3.30476146350799,-5.02349617178479,1.72712712683548,-5.7640062908429,-5.85636306964394,0.868398494423232,1.88830930227433,-6.2792515654936,-1.45446287725046,-3.44193680841572,1.40890069234139,1.97198384382082,-0.892799230269827,4.04321103907895,2.8004491138361,4.809418662581,-1.69278434137043,-3.07105250979535,-4.12728550731164,-3.59195392889046,4.84255143048937,-4.82819350420519,-1.53707134944857,-1.30423708029268,3.60712993650601,0.728050179732642,-4.80106774408418,-5.66750653981346,-1.509608509601,-5.51488674182709,0.776923059381931,-3.87263600516283,5.89158989063796,4.70331800564849,4.61113656050217,2.21314265946446,0.90042392406242,0.422243849305381,-2.70505383082815,5.30058887913081,2.20529708382031,-3.35191224868799,0.521812264060351,5.21624281598613,-1.5650143365568,-3.29878248230719,3.71086071593112,6.12350689261834,-4.27643013760671,-4.31705653647665,3.78649303516877,-1.3442081737974,1.12827342283103,-4.61236834920782,4.30649926544768,0.369446165323071,1.10545396096175,-4.27567585120588,-5.74200026407725,-6.04975287238822,-2.38380040821808,6.13615614816931,-3.6170130422302,0.767764068836605,1.7693249579597,-0.66509429706341,4.34615067039567,2.79919570905327,-2.10805145547397,-4.63789949353886,-3.03114074639526,-0.709760222300249,-5.92986375928032,-0.606308158022983,6.07891698060364,3.10456820820548,0.449921379333034,5.13057050308964,4.16303063432087,0.671851306742729,-6.22824668758768,-4.65364606845421,1.20361295950084,-3.56570177076684,3.84216147932214,-5.67021390374112,-3.57921212329167,-1.33556976556778,-3.662704554947,-3.03802718736891,-1.102137437956,0.236680242686595,-3.18505623919961,1.56403472566548,-5.27874089432847,3.99001304220007,-5.38424477724961,5.3505949760672,0.506023545401666,-1.20911103272357,-5.57048971747165,3.7580679990064,1.36593524963193,-5.21716827675198,-3.13142536486825,4.16004571972815,4.17058513145351,1.48754818516882,3.00743101696611,2.05043056577438,-4.30561948408442,1.3636475904741,3.67996969732016,3.18117867541642,5.77985295861158,0.464154835545818,3.7941499716753,-4.08254437053176,1.55975754396833,-6.1517397832717,-0.837386352048998,-1.19876123273315,0.368125659414894,2.2607426087514,-5.91791171291935,1.37256996508642,-0.0324295549002085,-5.01897129016896,0.439979735301953,-5.80959121535022,0.0551028771074674,1.1526752178303,4.23166198365619,-1.86376838376803,2.21869214107832,-5.18294858706374,-3.98690786989156,0.106092072531834,2.58778479810508,5.30370834707456,-4.12666256754546,4.56535051987233,2.08289267645298,5.8364923297747,1.46334398943707,1.57956042921289,0.792253094713655,-1.02569123888763,3.78020126453281,0.695629299936499,-0.894245822159327,-3.34037029336787,5.78005327305502,5.75706503725557,5.20355751538353,6.14532676013566,0.846449796490317;
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
