#ifndef __INIT_SHADERS__    
#define __INIT_SHADERS__ 1

#include <GL/glew.h>

#if defined (__APPLE__) || defined(MACOSX)
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif

// Define a helpful macro for handling offsets into buffer objects
#define BUFFER_OFFSET( offset )   ((GLvoid*) (offset))

typedef struct Shader   {   const char*  filename;
                            GLenum       type;
                            GLchar*      source;
                        }  tShader;

static char* readShaderSource(const char* shaderFile);

GLuint InitShader(const char* vShaderFile, const char* fShaderFile);

#endif //__INIT_SHADERS__   
