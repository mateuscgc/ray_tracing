uniform sampler2D tex;

//varying vec4 vColor;
varying vec2 vTexCoord;

void main() {
    // gl_FragColor = vColor;
    gl_FragColor = texture(tex, vTexCoord);
}
