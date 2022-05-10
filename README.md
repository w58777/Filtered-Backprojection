# Filtered-Backprojection

This project develops a CT reconstruction simulation software based on filtered back-projection algorithm, which can generate virtual X-rays for projection, and use R-L filter or S-L filter for filtering to generate reconstructed images.

- Program using C++ programming, using Qt for interface design and software release.
- The core steps in the program, such as projection, filtering and back-projection use self-programmed functions.
- The software can run in different computers without a specific environment.
- The program is extended from parallel beam to fan beam.

# Program description

The project file of the software is compiled based on qmake, written in Qt Creator, and managed by FBDemo.pro. The MainWindow class is defined in fbdemo.h, which contains the core functions of the software functions, such as projection, filtering-back-projection, etc. The implementation of the core functions is completed in fbdemo.cpp, and the fbdemo.ui file is used for interface design. 

You can load the project use Qt creator to modify.

# Usage

1. Open MyCTRS.exe, and directly to run
2. First click the Load Image button, select an image, support *.jpg, *.png, *.bmp format
3. Then select Parallel Projections or Fan Projections for projection;
4. Then select S-L filter or R-L filter for filtering;
5. Finally click Reconstruction to rebuild.

6. Precautions: Please follow the above steps. Due to limited time and insufficient program robustness, incorrect key sequence may cause program errors.

