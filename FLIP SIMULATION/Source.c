#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <GLFW/glfw3.h>

#define gX 0.0f
#define gY -9.81f
#define r 2.2f
#define h (r*2.2f)
#define particleNUM 200
#define gridX 100
#define gridY 100
#define cellX (int)(gridX/h)
#define cellY (int)(gridY/h)
#define dt (1.0f/60.0f)
#define damping 0.2f
#define k 0.7f //stiffness constant
#define repulsion .06f
#define alpha .9f //1 means pure flip, 0 means pure pic (flip -> pic)
#define overRelaxation 1.0f
#define rho0 1000
#define epsilon 0.000001

//hashing
#define SPATIAL_CELL_SIZE (2.2f * r)  // Slightly larger than particle diameter
#define SPATIAL_GRID_X ((int)(gridX / SPATIAL_CELL_SIZE) + 1)
#define SPATIAL_GRID_Y ((int)(gridY / SPATIAL_CELL_SIZE) + 1)
//hashing

//particle
float* particlePos = NULL;
float* particleVel = NULL;
//particle

//cell
int* cellType = NULL;
float* u = NULL;
float* v = NULL;
float* pu = NULL;
float* pv = NULL;   
float* du = NULL;
float* dv = NULL;
int* s = NULL;
float* divergence = NULL;
float* density = NULL;
float restDensity = 0.0f;
//cell

//spatial hash
int* spatialCellCount = NULL;
int* spatialCellStart = NULL;
int* spatialParticleIds = NULL;
//spatial hash

void spawn_particles() {
    int particlesPerRow = (int)sqrt(particleNUM);
    float space = 1.0f;

    float cubeWidth = particlesPerRow * space;
    float cubeHeight = ceil((float)particleNUM / particlesPerRow) * space;

    float startX = (gridX - cubeWidth) / 2.0f;
    float startY = (gridY - cubeHeight) / 2.0f;

    int index = 0;
    for (int y = 0; index < particleNUM; y++) {
        for (int x = 0; x < particlesPerRow && index < particleNUM; x++) {
            float px = startX + x * space;
            float py = startY + y * space;

            particlePos[index * 2 + 0] = px; // x
            particlePos[index * 2 + 1] = py; // y

            index++;
        }
    }
}

int cellCount = cellX*cellY;

void allocateMemory() {
// particles
particlePos = (float*)calloc(particleNUM * 2, sizeof(float)); // x,y
particleVel = (float*)calloc(particleNUM * 2, sizeof(float)); // vx,vy

// cells (Nx * Ny grid)
int numSpatialCells = SPATIAL_GRID_X * SPATIAL_GRID_Y;

cellType = (int*)calloc(cellCount, sizeof(int));
u = (float*)calloc(cellCount, sizeof(float));
v = (float*)calloc(cellCount, sizeof(float));
pu = (float*)calloc(cellCount, sizeof(float));
pv = (float*)calloc(cellCount, sizeof(float));
du = (float*)calloc(cellCount, sizeof(float));
dv = (float*)calloc(cellCount, sizeof(float));
s = (int*)calloc(cellCount, sizeof(float));
divergence = (float*)calloc(cellCount, sizeof(float)); // Updated variable name
density = calloc(cellCount, sizeof(float));
memset(s, 1, cellCount * sizeof(float));

spatialCellStart = (int*)calloc(numSpatialCells+1, sizeof(int));
spatialParticleIds = calloc(particleNUM, sizeof(int));
spatialCellCount = calloc(numSpatialCells, sizeof(int));
//spawnParticlesSquare(gridX * 0.5f, gridY * 0.5f, 40.0f);
}

void reset_Memory() {
    for (int i = 0; i < cellCount; i++) {
        density[i] = 0.0f;

    }
}
//void drawParticles() {
//    glPointSize(4.0f); // pixel size of particles
//    glBegin(GL_POINTS);
//
//    for (int i = 0; i < particleNUM; i++) {
//        float x = particlePos[2 * i];
//        float y = particlePos[2 * i + 1];
//
//        // Normalize to [-1,1] for OpenGL
//        float nx = (x / gridX) * 2.0f - 1.0f;
//        float ny = (y / gridY) * 2.0f - 1.0f;
//
//        glColor3f(0.2f, 0.6f, 1.0f); // blue-ish color
//        glVertex2f(nx, ny);
//    }
//
//    glEnd();
//}

void integrateParticles(int integrate) {
    for (int i = 0; i < particleNUM; i++) {
        // Apply gravity to velocity
        if (integrate) {
            particleVel[2 * i] += gX * dt;
            particleVel[2 * i + 1] += gY * dt;

            // Update positions
            particlePos[2 * i] += particleVel[2 * i] * dt;
            particlePos[2 * i + 1] += particleVel[2 * i + 1] * dt;

        }

        // Wall collisions
        //float x = particlePos[2 * i];
        //float y = particlePos[2 * i + 1];

        //// Right wall
        //if (x > gridX - r) {
        //    particlePos[2 * i] = gridX - r;
        //    particleVel[2 * i] *= -damping;
        //}
        //// Left wall
        //if (x < r) {
        //    particlePos[2 * i] = r;
        //    particleVel[2 * i] *= -damping;
        //}
        //// Top wall
        //if (y > gridY - r) {
        //    particlePos[2 * i + 1] = gridY - r;
        //    particleVel[2 * i + 1] *= -damping;
        //}
        //// Bottom wall
        //if (y < r) {
        //    particlePos[2 * i + 1] = r;
        //    particleVel[2 * i + 1] *= -damping;
        //}
    }
}



//delete later
void freeMemory() {
free(particlePos);
free(particleVel);

free(cellType);
free(u);
free(v);
free(pu);
free(pv);
free(divergence);
free(spatialCellCount);
free(spatialCellStart);
}
//delete later

float clamp(float x, float minVal, float maxVal) {
    if (x < minVal) return minVal;
    if (x > maxVal) return maxVal;
    return x;
}

void pushParticlesApart(int iter_) {
    float minDist = 2.0f * r;
    float minDist2 = minDist * minDist;

    int spatialGridX = SPATIAL_GRID_X;
    int spatialGridY = SPATIAL_GRID_Y;
    int numSpatialCells = spatialGridX * spatialGridY;

    for (int iter = 0; iter < iter_; iter++) {
        // Reset cell counts
        memset(spatialCellCount, 0, numSpatialCells * sizeof(int));

        // Count particles per cell
        for (int i = 0; i < particleNUM; i++) {
            float x = particlePos[2 * i];
            float y = particlePos[2 * i + 1];

            int xi = (int)(x / SPATIAL_CELL_SIZE);
            int yi = (int)(y / SPATIAL_CELL_SIZE);
            xi = clamp(xi, 0, spatialGridX - 1);
            yi = clamp(yi, 0, spatialGridY - 1);

            int cellIdx = xi * spatialGridY + yi;
            spatialCellCount[cellIdx]++;
        }


        // Build prefix sum
        //im using an inclusive bucket storage for the prefix sum
        int sum = 0;
        for (int i = 0; i < numSpatialCells; i++) {
            sum += spatialCellCount[i];
            spatialCellStart[i] = sum;
            //printf("sum: %d\n", spatialCellStart[i]);
        }
        spatialCellStart[numSpatialCells] = sum;

        // Reset counts for filling
        memset(spatialCellCount, 0, numSpatialCells * sizeof(int));

        // Assign particles to cells
        for (int i = 0; i < particleNUM; i++) {
            float x = particlePos[2 * i];
            float y = particlePos[2 * i + 1];

            int xi = (int)(x / SPATIAL_CELL_SIZE);
            int yi = (int)(y / SPATIAL_CELL_SIZE);
            xi = clamp(xi, 0, spatialGridX - 1);
            yi = clamp(yi, 0, spatialGridY - 1);

            int cellIdx = xi * spatialGridY + yi;
            int index = spatialCellStart[cellIdx] + spatialCellCount[cellIdx]++;
            spatialParticleIds[index] = i;
            //spatialCellCount[cellIdx]++;
        }

        // Resolve collisions
        for (int i = 0; i < particleNUM; i++) {
            float px = particlePos[2 * i];
            float py = particlePos[2 * i + 1];

            int pxi = (int)(px / SPATIAL_CELL_SIZE);
            int pyi = (int)(py / SPATIAL_CELL_SIZE);

            // Check 3x3 neighborhood
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    int xi = pxi + dx;
                    int yi = pyi + dy;

                    if (xi < 0 || xi >= spatialGridX || yi < 0 || yi >= spatialGridY) continue;

                    int cellIdx = xi * spatialGridY + yi;
                    int first = spatialCellStart[cellIdx];
                    int last = first + spatialCellCount[cellIdx];

                    for (int j = first; j < last; j++) {
                        int id = spatialParticleIds[j];
                        if (id == i) continue;

                        float qx = particlePos[2 * id];
                        float qy = particlePos[2 * id + 1];

                        float dx = qx - px;
                        float dy = qy - py;
                        float d2 = dx * dx + dy * dy;

                        if (d2 > minDist2 || d2 == 0.0f) continue;

                        float d = sqrtf(d2);
                        float s = repulsion * (minDist - d) / d;
                        dx *= s;
                        dy *= s;

                        particlePos[2 * i] -= dx;
                        particlePos[2 * i + 1] -= dy;
                        particlePos[2 * id] += dx;
                        particlePos[2 * id + 1] += dy;
                    }
                }
            }
        }
    }
}


//now we compute cell-particle density as rho
//lazy right now ill do transfer velocities and solve at a later date
void computeDensity() {
    for (int den = 0; den < cellCount; den++) {
        density[den] = 0.0f;
    }
    float h1 = 1.0f / h;
    float h2 = 0.5f * h;

    for (int i = 0; i < particleNUM; i++) {
        float x = clamp(particlePos[i * 2], h, (cellX-1)*h);
        float y = clamp(particlePos[i * 2 + 1], h, (cellY - 1) * h);

        int x0 = (int)((x - h2) * h1);
        float tx = ((x - h2) - x0 * h) * h1;
        int x1 = (int)min(x0 + 1, cellX - 2);

        int y0 = (int)((y - h2) * h1);
        float ty = ((y - h2) - y0 * h) * h1;
        int y1 = (int)min(y0 + 1, cellY - 2);

        float sx = 1.0f - tx;
        float sy = 1.0f - ty;

        if (x0 < cellX && y0 < cellY) density[x0 * cellY + y0] += sx * sy;
        if (x1 < cellX && y0 < cellY) density[x1 * cellY + y0] += tx * sy;
        if (x1 < cellX && y1 < cellY) density[x1 * cellY + y1] += tx * ty;
        if (x0 < cellX && y1 < cellY) density[x0 * cellY + y1] += sx * ty;
    }

    if (restDensity == 0.0f) {
        float sum = 0.0f;
        int numFluidCells = 0;
        int numCells = cellX * cellY;
        for (int cell = 0; cell < numCells; cell++) {
            if (cellType[cell] == 2) {
                sum += density[cell]; //if fluid compute density sum of cell;
                numFluidCells++;
            }
        }

        if (numFluidCells > 0) {
            restDensity = sum / numFluidCells;
        }
    }
}

void transferVelocity(int toGrid) {
    int ny = cellY;
    int nx = cellX;
    float h1 = 1.0f / h;
    float h2 = 0.5f * h;

    //reset cell
    if (toGrid) {
        memcpy(pu, u, sizeof(float) * cellCount);
        memcpy(pv, v, sizeof(float) * cellCount);
        for (int res = 0; res < cellCount; res++) {
            u[res] = 0.0f;
            v[res] = 0.0f;
            du[res] = 0.0f;
            dv[res] = 0.0f;
        }


        for (int i = 0; i < cellCount; i++) {
            cellType[i] = s[i] == 0 ? 0 : 1; //solid : air
        }

        for (int j = 0; j < particleNUM; j++) {
            float x = particlePos[j * 2];
            float y = particlePos[j * 2 + 1];
            int xi = (int)clamp(floor(x*h1),0.0f, nx - 1);
            int yi = (int)clamp(floor(y*h1),0.0f, ny - 1);
            int index = xi * ny + yi;
            if (cellType[index] == 1) cellType[index] = 2; // if air, make fluid type
        }
    }


    for (int comp = 0; comp < 2; comp++) {
        float dx = comp == 0 ? 0.0f : h2;
        float dy = comp == 0 ? h2 : 0.0f;

        float* f = comp == 0 ? u : v;
        float* prevF = comp == 0 ? pu : pv;
        float* d = comp == 0 ? du : dv;

        //now we do grid to particles
        //find 4 cells 
        for (int p = 0; p < particleNUM; p++) {
            float x = particlePos[p * 2];
            float y = particlePos[p * 2 + 1];

            x = clamp(x, h, (float)((nx - 1) * h));
            y = clamp(y, h, (float)((ny - 1) * h));

            int x0 = (int)clamp(floorf(x - dx), h, (float)((nx - 2)));
            int y0 = (int)clamp(floorf(y - dy), h, (float)((ny - 2)));
            //now we have cell coords

            //locate neighbor x
            //locate right and top cells
            int x1 = (int)min(x0 + 1, cellX - 2);
            int y1 = (int)min(y0 + 1, cellY - 2);

            //compensate stagger
            float tx = ((x - dx) - x0 * h) * h1;
            float ty = ((y - dy) - y0 * h) * h1;

            float sx = 1.0f - tx;
            float sy = 1.0f - ty;
            // compute weights

            float w0 = sx * sy;
            float w1 = tx * sy;
            float w2 = tx * ty;
            float w3 = sx * ty;

            int nr0 = x0 * cellY + y0;
            int nr1 = x1 * cellY + y0;
            int nr2 = x1 * cellY + y1;
            int nr3 = x0 * cellY + y1;


            if (toGrid) {
                float pv = particleVel[2 * p + comp];
                f[nr0] += pv * w0; d[nr0] += w0;
                f[nr1] += pv * w1; d[nr1] += w1;
                f[nr2] += pv * w2; d[nr2] += w2;
                f[nr3] += pv * w3; d[nr3] += w3;
            }
            else {
                // G2P transfer
                int offset = comp == 0 ? gridY : 1;
                float f0 = ((cellType[nr0] != 1) || cellType[nr0 - offset] != 1) ? 1.0f : 0.0f;
                float f1 = ((cellType[nr1] != 1) || cellType[nr1 - offset] != 1) ? 1.0f : 0.0f;
                float f2 = ((cellType[nr2] != 1) || cellType[nr2 - offset] != 1) ? 1.0f : 0.0f;
                float f3 = ((cellType[nr3] != 1) || cellType[nr3 - offset] != 1) ? 1.0f : 0.0f;
                float d = f0 * w0 + f1 * w1 + f2 * w2 + f3 * w3;
                float vel = particleVel[p * 2 + comp];


                // blend FLIP and PIC
                //particleVel[2 * p + comp] = (1.0f - alpha) * flip + alpha * pic;
                if (d > 0.0f) {
                    float pic = (f0 * w0 * f[nr0] + f1*w1*f[nr1] + f2 * w2 * f[nr2] + f3 * w3 * f[nr3])/d;
                    float corr = (
                        (f0 * w0 * (f[nr0] - prevF[nr0])) +
                        (f1 * w1 * (f[nr1] - prevF[nr1])) +
                        (f2 * w2 * (f[nr2] - prevF[nr2])) +
                        (f3 * w3 * (f[nr3] - prevF[nr3]))
                        )/d;
                    float flip = vel + corr;
                    particleVel[2 * p + comp] = alpha * flip + (1.0f - alpha) * pic;
                }
            }
        }
        if (toGrid) {
            for (int i = 0; i < cellCount; i++) {
                if (d[i] > 0.0f) {
                    f[i] /= d[i];
                }
            }
            for (int i = 0; i < cellX; i++) {
                for (int j = 0; j < cellY; j++) {
                    int solid = cellType[i * cellY + j];
                    if (solid || i > 0 && cellType[(i - 1) * cellY + j] == 0) {
                        u[i * cellY + j] = pu[i * cellY + j];
                        v[(i - 1) * cellY + j] = 0.0f;
                        u[(i - 1) * cellY + j] = 0.0f;
                    }

                    if (solid || j > 0 && cellType[i * cellY + j - 1] == 0) {
                        v[i * cellY + j] = pv[i * cellY + j];
                        v[i * cellY + j - 1] = 0.0f;
                        u[i * cellY + j - 1] = 0.0f;
                    }

                }
            }
        }
    }
}

void solveIncompressibility(int numIter) {
    memset(divergence, 0.0f, cellCount * sizeof(float));
    memcpy(pu, u, cellCount * sizeof(float));
    memcpy(pv, v, cellCount * sizeof(float));
    //reset divergence array and clone the previous velocity components for differences later
    float cp = rho0 * h / dt;
    //run based on user defined divergence/pressure solve iterations
    for (int iter = 0; iter < numIter; iter++) {
        for (int i = 1; i < cellX - 1; i++) {
            for (int j = 1; j < cellY - 1; j++) {
                if (cellType[i * cellY + j] == 0) continue;
                
                int center = i * cellY + j;
                int left = (i - 1) * cellY + j;
                int right = (i + 1) * cellY + j;
                int top = i * cellY + j + 1;
                int bottom = i * cellY + j - 1;
                //defined direct neighbors from center;

                int sc = s[center];
                int sl = s[left];
                int sr = s[right];
                int st = s[top];
                int sb = s[bottom];
                int sValidNum = sl + sr + st + sb;
                if (sValidNum == 0) continue;
                //validity

                //solve for divergence;
                float div = u[right] - u[center] + v[top] - v[center];
                
                if (restDensity > 0.0f) {
                    float compression = density[i * cellY + j] - restDensity;
                    if (compression > 0.0f) {
                        div -= k * compression;
                    }
                }

                float p = (-div / sValidNum)*overRelaxation;
                divergence[center] += cp * p;
                u[center] -= sl * p;
                u[right] += sr * p;
                v[top] += st * p;
                v[bottom] -= sb * p;

                
            }
        }
    }
}


//void solveIncompressibility(int numIter) {
//    float scale = dt / (rho0 * h * h);
//
//    for (int iter = 0; iter < numIter; iter++) {
//        for (int i = 1; i < cellX - 1; i++) {
//            for (int j = 1; j < cellY - 1; j++) {
//                int idx = i * cellY + j;
//                if (cellType[idx] != 2) continue; // Only fluid cells
//
//                int left = (i - 1) * cellY + j;
//                int right = (i + 1) * cellY + j;
//                int bottom = i * cellY + (j - 1);
//                int top = i * cellY + (j + 1);
//
//                // Count valid fluid neighbors
//                int validCount = 0;
//                if (cellType[left] == 2) validCount++;
//                if (cellType[right] == 2) validCount++;
//                if (cellType[top] == 2) validCount++;
//                if (cellType[bottom] == 2) validCount++;
//
//                if (validCount == 0) continue;
//
//                // Calculate divergence
//                float div = u[right] - u[idx] + v[top] - v[idx];
//
//                // Add density constraint - this is what makes it liquid-like
//                if (restDensity > 0.0f && density[idx] > 0.0f) {
//                    float densityError = (density[idx] - restDensity) / restDensity;
//                    div += 2.0f * densityError; // Strong density constraint
//                }
//
//                float pressure = -div / validCount;
//                pressure *= scale;
//
//                // Apply pressure gradient
//                u[idx] -= pressure;
//                u[right] += pressure;
//                v[idx] -= pressure;
//                v[top] += pressure;
//            }
//        }
//    }
//}

//void renderGrid() {
//    glColor3f(1.0f, 1.0f, 1.0f); // gray outlines
//    glLineWidth(1.0f);
//
//    for (int i = 0; i < gridX; i++) {
//        for (int j = 0; j < gridY; j++) {
//            float x0 = i * h;
//            float y0 = j * h;
//            float x1 = x0 + h;
//            float y1 = y0 + h;
//
//            // Convert to normalized OpenGL coordinates [-1,1]
//            float nx0 = (x0 / gridX) * 2.0f - 1.0f;
//            float ny0 = (y0 / gridY) * 2.0f - 1.0f;
//            float nx1 = (x1 / gridX) * 2.0f - 1.0f;
//            float ny1 = (y1 / gridY) * 2.0f - 1.0f;
//
//            glBegin(GL_LINE_LOOP);
//            glVertex2f(nx0, ny0);
//            glVertex2f(nx1, ny0);
//            glVertex2f(nx1, ny1);
//            glVertex2f(nx0, ny1);
//            glEnd();
//        }
//    }
//}

void setSolidCell(int i, int j) {
    int idx = i * cellY + j;
    cellType[idx] = 0;
    s[idx] = 0.0;  // make sure "validity" array says no fluid passes through
}

void setBoundaryWalls() {
    for (int i = 0; i < cellX; i++) {
        for (int j = 0; j < cellY; j++) {
            if (i == 0 || j == 0 || i == cellX - 1 || j == cellY - 1) {
                setSolidCell(i, j);
            }
        }
    }
}



int main() {
    allocateMemory();
    spawn_particles();
    //init + version declares
    GLFWwindow* window;
    if (!glfwInit()) {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return -1;
    }

    window = glfwCreateWindow(800, 800, "Fluid Sim", NULL, NULL);
    if (!window) {
        fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        return -1;
    }
    
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // vsync

    // --- OpenGL 2D setup ---
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    int count = 0;
    setBoundaryWalls();
    double lastTime = glfwGetTime();
    while (!glfwWindowShouldClose(window)) {
        //glfwPollEvents();
        //logic
        double ctime = glfwGetTime();
        double deltaTime = ctime - lastTime;
        lastTime = ctime;
        printf("frame: %d\ntime: %.2f\nframe/sec: %.2f\n", count++, ctime, count/ctime);
        integrateParticles(1);
        pushParticlesApart(5);
        integrateParticles(0);
        transferVelocity(1);
        computeDensity();
        solveIncompressibility(12);
        transferVelocity(0);

        //logic


        //boundary / collisions
        //boundary / collisions

       // --- Rendering ---
        glClear(GL_COLOR_BUFFER_BIT);
        //renderGrid();
        glLoadIdentity();
        //glColor3f(1.0f, 1.0f, 1.0f); // White particles
        glPointSize(3.5f);

        // In your rendering code
        glBegin(GL_POINTS);
        for (int i = 0; i < particleNUM; i++) {
            glColor3f(0.0f, 1.0f, 0.0f);

            float x = particlePos[i *2];
            float y = particlePos[i * 2 + 1];
            float nx = (x / gridX) * 2.0f - 1.0f;
            float ny = (y / gridY) * 2.0f - 1.0f;

            glVertex2f(nx, ny);
        }
        glEnd();
        glfwSwapBuffers(window);
        //printf("problem\n");
        glfwPollEvents();
        //drawParticles();
        //glfwSwapBuffers(window);
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;

}