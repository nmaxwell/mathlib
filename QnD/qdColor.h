#ifndef qdColor_H
#define qdColor_H




struct qdColor
{
public:
    uint8_t red;
    uint8_t green;
    uint8_t blue;
    
    qdColor(float R,float G,float B) // [0,1]_float -> [0,255]_int
    {
        if (R<0.f) R = 0.f; if (R>1.f) R = 1.f;
        if (G<0.f) G = 0.f; if (G>1.f) G = 1.f;
        if (B<0.f) B = 0.f; if (B>1.f) B = 1.f;    
        red = uint8_t(R*255);
        green = uint8_t(G*255);
        blue = uint8_t(B*255);
    }
    
    qdColor(uint8_t const & R, uint8_t const & G, uint8_t const & B) 
    :red(R),green(G),blue(B) {}
    
public:
    qdColor() {};
    ~qdColor() {};
    void operator= (const qdColor & rhs) {red = rhs.red; green = rhs.green; blue = rhs.blue;}
    qdColor(const qdColor & rhs) {operator=(rhs);}
};










#endif
