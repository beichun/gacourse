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

// Include GLM
#include <glm/glm.hpp>
GLFWwindow* window;

//define object shape
const int a = 1;
const int b = 1;
const int c = 1;

//define some parameters of cubes
const double Length = 0.1;

typedef Eigen::Array<double,Eigen::Dynamic,3,Eigen::RowMajor> Matrix3d;
typedef Eigen::Array<double,Eigen::Dynamic,2,Eigen::RowMajor> Matrix2d;

Matrix2d initBaseSpringtoMass() {
    Matrix2d baseSpringtoMass(28,2);
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

Matrix3d initBaseMassPosition(){
    Matrix3d baseMassPosition(8,3);
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

Matrix3d initMassPosition(){
    Matrix3d massPosition((a+1)*(b+1)*(c+1),3);
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
    //std::cout<<massPosition<<std::endl;
    return massPosition;
}

Matrix3d initMassVelocity(){
    Matrix3d massVelocity((a+1)*(b+1)*(c+1),3);
    for(int i=0;i<(a+1)*(b+1)*(c+1);i++){
        massVelocity.row(i) << 0,0,0;
    }
    //std::cout<<"Vel"<<std::endl<<massVelocity<<std::endl;
    return massVelocity;
}

Matrix3d initMassAcceleration(){
    Matrix3d massAcceleration((a+1)*(b+1)*(c+1),3);
    for(int i=0;i<(a+1)*(b+1)*(c+1);i++){
        massAcceleration.row(i) << 0,0,0;
    }
    //std::cout<<massAcceleration<<std::endl;
    return massAcceleration;
}

Matrix3d baseMassPosition = initBaseMassPosition();
Matrix3d massPosition = initMassPosition();
Matrix3d massVelocity = initMassVelocity();
Matrix3d massAcceleration = initMassAcceleration();
Matrix2d baseSpringtoMass = initBaseSpringtoMass();

Matrix2d initSpringtoMass(){
    Matrix2d springtoMass(28*a*b*c,2);
    for(int i=0;i<28*a*b*c;i++) {
        int cube_index = i / 28;
        //std::cout<<cube_index<<"cube_index"<<i<<std::endl;
        int spring_position = i % 28;
        int deltax = (cube_index % (a * b)) % a;
        int deltay = (cube_index % (a * b) / a);
        int deltaz = (cube_index / (a * b));
        //std::cout << "deltax=" << deltax << std::endl;
        //std::cout << "deltay=" << deltay << std::endl;
        //std::cout << "deltaz=" << deltaz << std::endl;
        for (int j = 0; j < 2; j++) {
            int mass_position = baseSpringtoMass(spring_position,j);
            int x0 = int(baseMassPosition(mass_position,0)/Length);
            int y0 = int(baseMassPosition(mass_position,1)/Length);
            int z0 = int(baseMassPosition(mass_position,2)/Length);
            springtoMass(i,j) = (z0 + deltaz) * (a + 1) * (b + 1) + (y0 + deltay) * (b + 1) + x0 + deltax;
        }
    }
    //std::cout << "springtoMass"<<std::endl<<springtoMass<<std::endl;
    return springtoMass;
}

Matrix2d springtoMass = initSpringtoMass();

int main() {

#pragma omp parallel
    printf("Hello, world.\n");

    // Initialise GLFW
    glewExperimental = true; // Needed for core profile
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL

    // Open a window and create its OpenGL context
    //GLFWwindow* window; // (In the accompanying source code, this variable is global for simplicity)
    window = glfwCreateWindow( 1024, 768, "Tutorial 01", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window); // Initialize GLEW
    glewExperimental=true; // Needed in core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return -1;
    }

    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

    // Dark blue background
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

    do{
        // Clear the screen. It's not mentioned before Tutorial 02, but it can cause flickering, so it's there nonetheless.
        glClear( GL_COLOR_BUFFER_BIT );

        // Draw nothing, see you in tutorial 2 !

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();

    } // Check if the ESC key was pressed or the window was closed
    while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           glfwWindowShouldClose(window) == 0 );
 
	return 0;
}