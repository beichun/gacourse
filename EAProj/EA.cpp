// Include standard headers
#include <stdio.h>
#include <stdlib.h>

//include header for io
#include <iostream>

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


typedef Eigen::Array<double,Eigen::Dynamic,3,Eigen::RowMajor> ArrayX3dRowMajor;
typedef Eigen::Array<double,Eigen::Dynamic,2,Eigen::RowMajor> ArrayX2dRowMajor;

GLFWwindow* window;

//define object shape
const int a = 1;
const int b = 5;
const int c = 1;

const int num_masses = (a+1)*(b+1)*(c+1);
const int num_cubes = a*b*c;
const int num_springs = 28*num_cubes;

//define some parameters of cubes
const double PI = 3.1415926;
const double Mass = 0.1;  //kg
const double SpringConstraint = 1000;  //N/m
const double Length = 0.1;
const double w = 4*PI;
const double CoefaRange = 0.02;
const double CoefbRange = PI;

//friction coefficients for glass on glass
const double Muk = 0.8;
const double Mus = 1;

const double InitialHeight = 0.5;
const double InitialVelocityY = 0.0;
const double TimeStep = 0.0005;  //s
const double Gravity[3] = {0,0,-9.81};
//
//const double Gravity[3] = {0,0,-0.5};

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


Eigen::ArrayXd initSpringCoefa(){

    Eigen::ArrayXd springCoefa(num_cubes);
    springCoefa = CoefaRange*Eigen::ArrayXd::Random(num_cubes);
    return springCoefa;
}


Eigen::ArrayXd initSpringCoefb(){
    Eigen::ArrayXd springCoefb(num_cubes);
    springCoefb = CoefbRange*Eigen::ArrayXd::Random(num_cubes);
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
        Eigen::ArrayXd& springCoefa,
        Eigen::ArrayXd& springCoefb,
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
    window = glfwCreateWindow( 1024, 768, "Tutorial 04 - Colored Cube", NULL, NULL);
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
            glm::vec3(1,1.5,1.5),
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
        Eigen::ArrayXd& springCoefa,
        Eigen::ArrayXd& springCoefb,
        Eigen::ArrayXd& l0,
        ArrayX3dRowMajor& massPosition,
        ArrayX2dRowMajor& springtoMass,
        ArrayX3dRowMajor& massVelocity,
        ArrayX3dRowMajor& massAcceleration,
        ArrayX3dRowMajor& position_history){

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
    for (int i=0;i<c;i++){
        for (int j=(a+1)*b;j<(a+1)*(b+1);j++)
        fitness += massPosition(i*(a+1)*(b+1)+j,1);
    }
    fitness /= (a+1)*(c+1);
    //std::cout<<"fitness="<<fitness<<std::endl;
    return fitness;
}

void mutate(
        Eigen::ArrayXd& springCoefa_new,
        Eigen::ArrayXd& springCoefb_new){
    double r = ((double) rand() / (RAND_MAX));
    int r1 = int(2*r*num_cubes);
    if (r1<num_cubes) {
        springCoefa_new(r1) = CoefaRange*((double) rand() / (RAND_MAX));
    }
    else{
        springCoefb_new(r1-num_cubes) = CoefbRange*((double) rand() / (RAND_MAX));
    }
    //std::cout<<r1<<std::endl;
}

int main() {

    Eigen::ArrayXi cubeVertex = initCubeVertex();
    Eigen::ArrayXi triangleVertex = initTriangleVertex();
    ArrayX3dRowMajor baseMassPosition = initBaseMassPosition();
    ArrayX3dRowMajor massVelocity = initMassVelocity();
    ArrayX3dRowMajor massAcceleration = ArrayX3dRowMajor::Zero(num_masses,3);
    ArrayX2dRowMajor baseSpringtoMass = initBaseSpringtoMass();
    ArrayX2dRowMajor springtoMass = initSpringtoMass(
            baseSpringtoMass,
            baseMassPosition);
    Eigen::ArrayXd l0 = initL0();


    int num_frames = 5000;
    int skip_frames = 32;
    int num_evaluations = 1000;
    /*#pragma omp parallel
    printf("Hello, world.\n");*/

    double new_fitness;
    double best_fitness = 0;
    ArrayX3dRowMajor position_history = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);
    ArrayX3dRowMajor best_position_history = ArrayX3dRowMajor::Zero(num_frames*num_masses,3);

    Eigen::ArrayXd springCoefa_best = initSpringCoefa();
    Eigen::ArrayXd springCoefb_best = initSpringCoefb();

    Eigen::ArrayXd springCoefa_new = springCoefa_best;
    Eigen::ArrayXd springCoefb_new = springCoefb_best;

//    #pragma omp parallel for
    for (int i=0;i<num_evaluations;i++){
        // initialize the structure
        springCoefa_new = springCoefa_best;
        springCoefb_new = springCoefb_best;

        //mutate
        mutate(springCoefa_new,springCoefb_new);

        // simulate it get the fitness
        ArrayX3dRowMajor massPosition = initMassPosition();

        new_fitness = simulate(
                num_frames,
                skip_frames,
                springCoefa_new,
                springCoefb_new,
                l0,
                massPosition,
                springtoMass,
                massVelocity,
                massAcceleration,
                position_history);
        //std::cout<<new_fitness<<std::endl;
        if (new_fitness>=best_fitness) {
            best_fitness = new_fitness;
            best_position_history = position_history;
            springCoefa_best = springCoefa_new;
            springCoefb_best = springCoefb_new;
        }
        std::cout<<"best_fitness "<<i<<" = "<<best_fitness<<std::endl;
    }
    //std::cout<<"best_fitness "<<best_fitness<<std::endl;
    std::cout<<"springCoefa_best = \n"<<springCoefa_best<<std::endl;
    std::cout<<"springCoefa_best = \n"<<springCoefb_best<<std::endl;

    // render this in glfw
    render(best_position_history,
           triangleVertex,
           cubeVertex,
           springtoMass);


    return 0;
}