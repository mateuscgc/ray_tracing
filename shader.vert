attribute vec3 aPosition;
// attribute vec4 aColor;
attribute vec2 aTexCoord;

uniform mat4 uProjectionMatrix;
uniform mat4 uViewMatrix;
uniform mat4 uModelMatrix;

// varying vec4 vColor;
varying vec2 vTexCoord;

void main() {
    vTexCoord = aTexCoord;
    // vColor = aColor;
    gl_Position = uProjectionMatrix * uViewMatrix * uModelMatrix * vec4(aPosition, 1.0);
    
}
