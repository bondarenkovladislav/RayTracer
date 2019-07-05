#include <math.h>

#include <GLUT/glsmapint.h>
#include <GLUT/GLUT.h>//Необходимо изменить данные файлы библиотек в зависимости от вашей ОС
#include "Bondarenko.hpp"
#include <stdio.h>


using namespace std;

double ScaleX = 0;
double ScaleY = 0;
double ScaleZ = 0;

point hitPoint;
vecteur e1, e2, trace, lightTrace;
ray n;

double *x = new double[3];
int * gauss();
int *z = new int[3];
int count = 0;


void init (void) {glClearColor(0.0,0.0,0.0,0.0);}

int windowWidth = 500, windowHeight = 500;








//-------------------------Раздел параметров ------------------------------------

//Задание координат плоскости
vecteur e3={0.0,1.0,1.0};//Нормальный вектор
plane Plane {
    e3,
    //{100.0, 100.0, 50.0}//d
    {1,1,0}//Точка принадлежащая плоскости
};
//Задание координат сфер
const int cntSpheres = 2;
sphere Spheres[cntSpheres] =
{
    {
        {-100.0, 290.0, 0.0},//Центр
        100.0,//Размер
        0//id Материала
    },
    {
        {100.0,100.0,100.0},
        100,
        1
    }
};
//Задание материалов сфер
const int cntMaterials = 2;
const material Materials[cntMaterials] =
{
    {
        1,//Коэффицент отражения
        1.0, 0.0, 1.0//Цвет поверхности
    },
    //Аналогично для второго материала
    {
        0.5,
        1.0, 0.0, 0.0
    }
};
//Задание источника света
const int cntLights = 1;
light Lights[cntLights] =
{
    {
        {0.0, 0.0, -100.0},//pos
        1.0, 1.0, 1.0 //rgb
    }
};
// --------------------------------------------------------------------







bool isBasis(){
    double a11 = e1.x;
    double a12 = e2.x;
    double a13 = e3.x;
    double a21 = e1.y;
    double a22 = e2.y;
    double a23 = e3.y;
    double a31 = e1.z;
    double a32 = e2.z;
    double a33 = e3.z;
    double res = (a11*a22*a33)+(a12*a23*a31)+(a13*a21*a32)-(a13*a22*a31)-(a11*a32*a23)-(a12*a21*a33);
    if(res == 0)
        return false;
    return true;
}

void createBasis(plane Plane)
{
    e1 = {Plane.d.x * 2, Plane.d.y*2, Plane.d.z *2};
    e2 = {e1.y, -e1.x, 0};
}



/*
 const int cntLights = 2;
 const light Lights[cntLights] =
 {
 {
 {0.0, 240.0, -100.0},
 1.0, 1.0, 1.0
 },
 {
 {640.0, 240.0, -10000.0},
 0.6, 0.7, 1.0
 }
 };*/

bool hitPlane(plane& plane,ray& ray)
{
    float alpha = plane.normal * ray.dir;
    float result;
    
    
    if(alpha != 0.0f)
    {
        result = - ((plane.normal.x* ray.start.x + plane.d.x) + (plane.normal.y + ray.start.y + plane.d.y) + (plane.normal.z + ray.start.z + plane.d.z));
        result = result / alpha;
        hitPoint.x = ray.start.x + ray.dir.x * result;
        hitPoint.y = ray.start.y + ray.dir.y * result;
        hitPoint.z = ray.start.z + ray.dir.z * result;
        return (result >= 0.0f);
    }
    
    
    return false;
}

bool hitSphere(const ray &r, const sphere &s, float &t)
{
    vecteur dist = s.pos - r.start;
    float B = r.dir * dist;
    float D = B*B - dist * dist + s.size * s.size;
    if (D < 0.0f)
        return false;
    float t0 = B - sqrtf(D);
    float t1 = B + sqrtf(D);
    bool retvalue = false;
    if ((t0 > 0.1f) && (t0 < t))
    {
        t = t0;
        retvalue = true;
    }
    if ((t1 > 0.1f) && (t1 < t))
    {
        t = t1;
        retvalue = true;
    }
    return retvalue;
}

float degrees (vecteur v, vecteur n)
{
    float cos = (v * n)/(sqrtf(v.x*v.x + v.y * v.y +v.z * v.z));
    return (cos);
}

int * gauss()
{
    
    double a[3][3] = {e1.x,e1.y,e1.z,e2.x,e2.y,e2.z,e3.x,e3.y,e3.z};
    double y[3] = {lightTrace.x, lightTrace.y, lightTrace.z};
    if(!isBasis())
    {
        printf("Oh No");
        return NULL;
    }
    int *res = new int[3];
    double *x, max;
    int k, index;
    const double eps = 0.00001;  // точность
    x = new double[3];
    k = 0;
    while (k < 3)
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < 3; i++)
        {
            if (abs(a[i][k]) > max)
            {
                max = abs(a[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            return 0;
        }
        for (int j = 0; j < 3; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < 3; i++)
        {
            double temp = a[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < 3; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < 3; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = 2; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    
    x[0] = round(x[0]);
    x[1] = round(x[1]);
    x[2] = round(x[2]);
    
    res[0] = x[0];
    res[1] = x[1];
    res[2] = x[2];
    
    return res;
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    Spheres[0].pos.x += ScaleX;
    Spheres[0].pos.y += ScaleY;
    Spheres[0].pos.z += ScaleZ;
    glPushMatrix();
    glTranslatef(0, 0, -1000);
    
    for (int y = -500; y < windowHeight; ++y)
    {
        for (int x =-500; x < windowWidth; ++x)
        {
            float red = 0, green = 0, blue = 0;
            float coef = 1.0f;
            int level = 0;
            
            ray viewRay = { {float(x), float(y), -1000.0f}, { 0.0f, 0.0f, 1.0f}};
            do
            {
                float t = 2000.0f;
                int currentSphere= -1;
                for (int i = 0; i < cntSpheres; ++i)
                {
                    if (hitSphere(viewRay, Spheres[i], t))
                    {
                        currentSphere = i;
                    }
                }
                if (currentSphere == -1)
                {
                    if(level !=0 && (e3.y == 0 && e3.x == 0))
                        break;
                    int currentPlane= (hitPlane(Plane,viewRay) ? 1:-1);
                    if (currentPlane == -1) break;
                    trace  = Lights[0].pos - hitPoint;
                    n.start = {trace.x, trace.y,trace.z};
                    n.dir ={0.0,0.0,1.0};
                    hitPlane(Plane, n);
                    vecteur light = hitPoint - viewRay.start;
                    lightTrace = hitPoint - Plane.d;
                    createBasis(Plane);
                    z = gauss();
                    float temp =degrees(light, Plane.normal);
                    if(((z[0]+z[1])/50)%2 == 0)
                    //if((z[0]/15)%2 == 0 && (z[1]/15)%2 == 0 && (z[2]/15)%2 == 0)
                    {
                        //red= green = blue = 1; //pow(temp, 1);
                        red -= pow(temp, 1);
                        blue -= pow(temp, 1);
                        green -= pow(temp,1);
                        if (level == 0)
                        {
                            red = green = blue = pow(temp, 5);
                        }
                        count++;
                    }
                    else if(level != 0)
                    {
                        red += pow(temp,2);
                        green += pow(temp,2);
                        blue +=pow(temp,2);
                    }
                    break;
                }
                point newStart = viewRay.start + t * viewRay.dir;
                vecteur n = newStart - Spheres[currentSphere].pos;
                float temp = n * n;
                if (temp == 0.0f) break;
                
                temp = 1.0f / sqrtf(temp);
                n = temp * n;
                
                material currentMat = Materials[Spheres[currentSphere].materialId];
                for (int j = 0; j < cntLights; ++j)
                {
                    light current = Lights[j];
                    vecteur dist = current.pos - newStart;
                    if (n * dist <= 0.0f) continue;
                    
                    float t = sqrtf(dist * dist);
                    if ( t <= 0.0f ) continue;
                    
                    ray lightRay;
                    lightRay.start = newStart;
                    lightRay.dir = (1/t) * dist;
                    
                    bool inShadow = false;
                    for (int i = 0; i < cntSpheres; ++i)
                    {
                        if (hitSphere(lightRay, Spheres[i], t))
                        {
                            inShadow = true; break;
                        }
                    }
                    
                    if (!inShadow)
                    {
                        // lambert
                        float lambert = (lightRay.dir * n) * coef;
                        red += lambert * current.red * currentMat.red;
                        green += lambert * current.green * currentMat.green;
                        blue += lambert * current.blue * currentMat.blue;
                    }
                }
                //Просчет отразившегося луча
                coef *= currentMat.reflection;
                float reflet = 2.0f * (viewRay.dir * n);
                viewRay.start = newStart;
                viewRay.dir = viewRay.dir - reflet * n;
                level++;
            }
            while (level < 2);
            //while ((coef > 0.0f) && (level < 10));

            glPointSize(1.0);
            glBegin(GL_POINTS);
            glColor3f(red, green, blue);
            glVertex3f(x, y, 0);
            glEnd();
        }
    }
    glPopMatrix();
    ScaleY = 0;
    ScaleX = 0;
    ScaleZ = 0;
    glFlush();
}

void reshape(int w, int h)
{
    windowWidth = w; windowHeight = h;
    glViewport(0,0,(GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //gluOrtho2D(-windowWidth, windowWidth, -windowHeight, windowHeight);
    //gluPerspective(45, 1, 1,400);
    gluPerspective(45.0, 1, 0.1, 3000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void SKeyboard(int key, int x, int y)
{
    switch(key)
    {
        case GLUT_KEY_LEFT: ScaleX-=50;
            break;
        case GLUT_KEY_RIGHT: ScaleX+=50;
            break;
        case GLUT_KEY_UP: ScaleY+=50;
            break;
        case GLUT_KEY_DOWN: ScaleY-=50;
    }
}

void Timer(int value)
{
    glutPostRedisplay();
    glutTimerFunc(50, Timer, 0);
}

void Keyboard(unsigned char key, int x, int y)
{
    switch(key)
    {
        case '+': ScaleZ += 50;
            break;
        case '-': ScaleZ -= 50;
            break;
    }
}

int main(int argc, char **argv)
{
    printf("C помощью стрелок осуществляется движение сферы. С помощью кнопок '+' '-' осуществляется приближение и отдаление сферы.\n \nКоординаты центра сфер, луча и плоскости задаются выше. Для того, чтобы сферы не отражали плоскость, измените координаты нормального вектора плоскости на {0,0,1}.\n");
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB|GLUT_DEPTH);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(100,100);
    glutCreateWindow("Bondarenko");
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(50, Timer, 0);
    glutSpecialFunc(SKeyboard);
    glutKeyboardFunc(Keyboard);
    glutMainLoop();
    return 0;
}

