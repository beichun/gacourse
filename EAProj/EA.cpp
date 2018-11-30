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

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

GLFWwindow* window;

//define object shape
const int a = 1;
const int b = 1;
const int c = 1;

//define some parameters of cubes
const double Length = 0.1;

typedef Eigen::Array<double,Eigen::Dynamic,3,Eigen::RowMajor> ArrayX3dRowMajor;
typedef Eigen::Array<double,Eigen::Dynamic,2,Eigen::RowMajor> ArrayX2dRowMajor;

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
    ArrayX3dRowMajor massPosition((a+1)*(b+1)*(c+1),3);
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

ArrayX3dRowMajor initMassVelocity(){
    ArrayX3dRowMajor massVelocity((a+1)*(b+1)*(c+1),3);
    for(int i=0;i<(a+1)*(b+1)*(c+1);i++){
        massVelocity.row(i) << 0,0,0;
    }
    //std::cout<<"Vel"<<std::endl<<massVelocity<<std::endl;
    return massVelocity;
}

ArrayX3dRowMajor initMassAcceleration(){
    ArrayX3dRowMajor massAcceleration((a+1)*(b+1)*(c+1),3);
    for(int i=0;i<(a+1)*(b+1)*(c+1);i++){
        massAcceleration.row(i) << 0,0,0;
    }
    //std::cout<<massAcceleration<<std::endl;
    return massAcceleration;
}

ArrayX3dRowMajor baseMassPosition = initBaseMassPosition();
ArrayX3dRowMajor massPosition = initMassPosition();
ArrayX3dRowMajor massVelocity = initMassVelocity();
ArrayX3dRowMajor massAcceleration = initMassAcceleration();
ArrayX2dRowMajor baseSpringtoMass = initBaseSpringtoMass();

ArrayX2dRowMajor initSpringtoMass(){
    ArrayX2dRowMajor springtoMass(28*a*b*c,2);
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

ArrayX2dRowMajor springtoMass = initSpringtoMass();

int render(){
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
    glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

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
            glm::vec3(4,3,-3), // Camera is at (4,3,-3), in World Space
            glm::vec3(0,0,0), // and looks at the origin
            glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
    );
    // Model matrix : an identity matrix (model will be at the origin)
    glm::mat4 Model      = glm::mat4(1.0f);
    // Our ModelViewProjection : multiplication of our 3 matrices
    glm::mat4 MVP        = Projection * View * Model; // Remember, matrix multiplication is the other way around

    // Our vertices. Tree consecutive floats give a 3D vertex; Three consecutive vertices give a triangle.
    // A cube has 6 faces with 2 triangles each, so this makes 6*2=12 triangles, and 12*3 vertices
    static const GLfloat g_vertex_buffer_data[] = {
            -1.0f,-1.0f,-1.0f,
            -1.0f,-1.0f, 1.0f,
            -1.0f, 1.0f, 1.0f,
            1.0f, 1.0f,-1.0f,
            -1.0f,-1.0f,-1.0f,
            -1.0f, 1.0f,-1.0f,
            1.0f,-1.0f, 1.0f,
            -1.0f,-1.0f,-1.0f,
            1.0f,-1.0f,-1.0f,
            1.0f, 1.0f,-1.0f,
            1.0f,-1.0f,-1.0f,
            -1.0f,-1.0f,-1.0f,
            -1.0f,-1.0f,-1.0f,
            -1.0f, 1.0f, 1.0f,
            -1.0f, 1.0f,-1.0f,
            1.0f,-1.0f, 1.0f,
            -1.0f,-1.0f, 1.0f,
            -1.0f,-1.0f,-1.0f,
            -1.0f, 1.0f, 1.0f,
            -1.0f,-1.0f, 1.0f,
            1.0f,-1.0f, 1.0f,
            1.0f, 1.0f, 1.0f,
            1.0f,-1.0f,-1.0f,
            1.0f, 1.0f,-1.0f,
            1.0f,-1.0f,-1.0f,
            1.0f, 1.0f, 1.0f,
            1.0f,-1.0f, 1.0f,
            1.0f, 1.0f, 1.0f,
            1.0f, 1.0f,-1.0f,
            -1.0f, 1.0f,-1.0f,
            1.0f, 1.0f, 1.0f,
            -1.0f, 1.0f,-1.0f,
            -1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f,
            -1.0f, 1.0f, 1.0f,
            1.0f,-1.0f, 1.0f
    };

    // One color for each vertex. They were generated randomly.
    static const GLfloat g_color_buffer_data[] = {
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

    GLuint vertexbuffer;
    glGenBuffers(1, &vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

    GLuint colorbuffer;
    glGenBuffers(1, &colorbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_color_buffer_data), g_color_buffer_data, GL_STATIC_DRAW);

    do{

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
        glDrawArrays(GL_TRIANGLES, 0, 12*3); // 12*3 indices starting at 0 -> 12 triangles

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();

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

int main() {

    #pragma omp parallel
    printf("Hello, world.\n");
    // initialize the structure

    // simulate it get the fitness

    // render this in glfw
    render();




	return 0;
}