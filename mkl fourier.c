#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <GL/glut.h>
#include <math.h>

// 定义链表节点结构
typedef struct Node {
    float x, y;
    struct Node* next;
} Node;

// 定义链表结构，包含头尾指针
typedef struct {
    Node* head;
    Node* tail;
} LinkedList;

// 全局链表定义
LinkedList list1, list2, list3, list4;
int num = 2;
int dataSize = 0;
double a;
int b;
char* resultFilePath1;
char* resultFilePath2;

// 函数声明
Node* createNode(float x, float y);
void addNode(LinkedList* list, float x, float y);
void readdata(const char* filename);
void display(void);
void initwin(double a);
void fftAndDisplay();
void saveFFTResults(const char* filename, float* data, int size);

// 创建新节点
Node* createNode(float x, float y) {
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->x = x;
    newNode->y = y;
    newNode->next = NULL;
    return newNode;
}

// 添加节点到链表
void addNode(LinkedList* list, float x, float y) {
    Node* newNode = createNode(x, y);
    if (list->head == NULL) {
        list->head = newNode;
        list->tail = newNode;
    } else {
        list->tail->next = newNode;
        list->tail = newNode;
    }
}

// 读取数据
void readdata(const char* filename) {
    list1.head = NULL;
    list1.tail = NULL;
    list2.head = NULL;
    list2.tail = NULL;

    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("无法打开文件\n");
        return;
    }

    char buffer[1024];
    int row = 0;
    while (fgets(buffer, 1024, file)) {
        if (row > 0) {
            float x, ch1, ch2;
            sscanf(buffer, "%f,%f,%f", &x, &ch1, &ch2);
            addNode(&list1, x, ch1);
            addNode(&list2, x, ch2);
            num++;
        }
        row++;
    }
    dataSize = row - 1;
    fclose(file);
}

// 初始化OpenGL窗口
void initwin(double a) {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-0.08 * num, num/2.0, -0.1 * a, a);
}

// 绘图函数
void display(void) {
    glClear(GL_COLOR_BUFFER_BIT);

    // 绘制浅灰色坐标轴和网格线
    glColor3f(0.8, 0.8, 0.8);

    glBegin(GL_LINES);
    glVertex2f(0, -0.1 * a);
    glVertex2f(0, a);
    glEnd();

    glBegin(GL_LINES);
    glVertex2f(-0.1 * num, 0);
    glVertex2f(num, 0);
    glEnd();

    for (int j = 0; j < 10; j++) {
        float x = j * (num / 10.0);
        char label[10];
        sprintf(label, "%.1f", x);
        glRasterPos2f(x, -0.05 * a);
        for (char* c = label; *c != '\0'; c++) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
        }
        glBegin(GL_LINES);
        glVertex2f(x, 0);
        glVertex2f(x, -0.03 * a);
        glEnd();
    }

    for (int j = 0; j < a; j++) {
        glBegin(GL_LINES);
        glVertex2f(-0.03 * num, j);
        glVertex2f(0, j);
        glEnd();
        for (int i = 2; i < 10; i++) {
            glBegin(GL_LINES);
            glVertex2f(-0.02 * num, j + i / 10.0);
            glVertex2f(0, j + i / 10.0);
            glEnd();
        }
    }
    fftAndDisplay();
}

// 保存 FFT 结果到文件
void saveFFTResults(const char* filename, float* data, int size) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("无法打开文件以保存结果\n");
        return;
    }

    for (int i = 0; i < size; i++) {
        fprintf(file, "%f\n", fabs(data[i]));
    }

    fclose(file);
}

// 执行 FFT 并显示结果
void fftAndDisplay() {
    MKL_LONG n = dataSize;
    MKL_LONG status;
    DFTI_DESCRIPTOR_HANDLE hand;
    float* x1 = (float*)malloc(n * sizeof(float));
    float* x2 = (float*)malloc(n * sizeof(float));

    // 复制数据到 MKL 数组
    Node* temp1 = list1.head;
    Node* temp2 = list2.head;
    for (int i = 0; temp1 != NULL && temp2 != NULL && i < dataSize; i++) {
        x1[i] = temp1->y;
        x2[i] = temp2->y;
        temp1 = temp1->next;
        temp2 = temp2->next;
    }

    // 创建并执行 FFT 描述符
    status = DftiCreateDescriptor(&hand, DFTI_SINGLE, DFTI_REAL, 1, n);
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_INPLACE);
    status = DftiCommitDescriptor(hand);

    status = DftiComputeForward(hand, x2);
    saveFFTResults(resultFilePath2, x2, n);

    status = DftiComputeForward(hand, x1);
    saveFFTResults(resultFilePath1, x1, n);

    glPointSize(5.0);
    glColor3f(0.0, 0.0, 0.0); // 黑色
    glBegin(GL_POINTS);
    for (int i = 1; i < n/2 ;i++) {
        glVertex2f(i, fabs(sqrt(x2[i]*x2[i]+x2[i+1]*x2[i+1])));
    }
    glEnd();
    glPointSize(5.0);
    glColor3f(1.0, 0.0, 0.0); // 红色
    glBegin(GL_POINTS);
    for (int i = 1; i < n/2 ;i++) {
        glVertex2f(i, fabs(sqrt(x1[i]*x1[i]+x1[i+1]*x1[i+1])));
     }
    glEnd();

    glFlush();

    // 清理
    DftiFreeDescriptor(&hand);
    free(x1);
    free(x2);
}

int main(int argc, char** argv) {
    if (argc < 4) {
        printf("Usage: %s <datafile> <resultfile1> <resultfile2>\n", argv[0]);
        return -1;
    }

    readdata(argv[1]);
    resultFilePath1 = argv[2];
    resultFilePath2 = argv[3];
    a = 60.0;
    b = 1;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(800, 600);
    glutCreateWindow("FFT Visualization");
    initwin(a);
    glutDisplayFunc(display);
    glutMainLoop();

    return 0;
}

