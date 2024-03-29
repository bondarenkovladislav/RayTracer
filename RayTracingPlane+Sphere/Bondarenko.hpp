#ifndef __ARITH_H__
#define __ARITH_H__

const double PIOVER180 = 0.017453292519943295769236907684886;

struct point {
    float x, y, z;
};

struct vecteur {
    float x, y, z;
    
    vecteur& operator += (const vecteur &v2){
        this->x += v2.x;
        this->y += v2.y;
        this->z += v2.z;
        return *this;
    }
};

inline point operator + (const point&p, const vecteur &v){ //Встраиваемые функции
    point p2={p.x + v.x, p.y + v.y, p.z + v.z };
    return p2;
}

inline point operator - (const point&p, const vecteur &v){
    point p2={p.x - v.x, p.y - v.y, p.z - v.z };
    return p2;
}

inline vecteur operator + (const vecteur&v1, const vecteur &v2){
    vecteur v={v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
    return v;
}

inline vecteur operator - (const point&p1, const point &p2){
    vecteur v={p1.x - p2.x, p1.y - p2.y, p1.z - p2.z };
    return v;
}

inline vecteur operator * (float c, const vecteur &v)
{
    vecteur v2={v.x *c, v.y * c, v.z * c };
    return v2;
}

inline vecteur operator - (const vecteur&v1, const vecteur &v2){
    vecteur v={v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
    return v;
}

inline float operator * (const vecteur&v1, const vecteur &v2 ) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}



struct color {
    enum OFFSET //Перечисления
    {
        OFFSET_RED = 0,
        OFFSET_GREEN = 1,
        OFFSET_BLUE = 2,
        OFFSET_MAX  = 3
    };
    float red, green, blue;
    
    inline color & operator += (const color &c2 ) {
        this->red +=  c2.red;
        this->green += c2.green;
        this->blue += c2.blue;
        return *this;
    }
    
    inline float & getChannel(OFFSET offset )
    {
        return reinterpret_cast<float*>(this)[offset]; //Приведение типов
    }
    
    inline float getChannel(OFFSET offset ) const
    {
        return reinterpret_cast<const float*>(this)[offset];
    }
};

inline color operator * (const color&c1, const color &c2 ) {
    color c = {c1.red * c2.red, c1.green * c2.green, c1.blue * c2.blue};
    return c;
}

inline color operator + (const color&c1, const color &c2 ) {
    color c = {c1.red + c2.red, c1.green + c2.green, c1.blue + c2.blue};
    return c;
}

inline color operator * (float coef, const color &c ) {
    color c2 = {c.red * coef, c.green * coef, c.blue * coef};
    return c2;
}

struct material {
    float reflection;
    float red, green, blue;
};

struct sphere {
    point pos;
    float size;
    int materialId;
};

struct light {
    point pos;
    float red, green, blue;
};

struct ray {
    point start;
    vecteur dir;
};

struct plane{
    vecteur normal;
    point d;
};

#endif

