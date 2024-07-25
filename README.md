## Repository Structure
* Additional Part/ is the directory with content for the Additional Requirement. Submodules are:
  - Additional_version.cpp, parallel code with additional implement.
  - Test.cpp, one more additional implement that can be used to check the parallel result with original result.
* animation/ is the directory with the visualization of parallel result.
   - neumann_C1.gif, visualization of C1 with neumann boundary condition.
   - neumann_C2.gif, visualization of C2 with neumann boundary condition.
   - periodic_C1.gif, visualization of C1 with periodic boundary condition.
   - periodic_C2.gif, visualization of C2 with periodic boundary condition.
  * Serial_Reaction_Equations.cpp, the original code.
* benchmark.pbs, can be used to submit in hpc to run the original serial code and counting the execute time.
* evaluate.pbs, can be used to submit in hpc to run the parallel code and counting the execute time.
* parallel_version.cpp, the parallel version of the original code.
* visualization.py, the post-processing code, used to visualize and generate the animation using the parallel result.
* reference.txt, the reference list.
* Documentation.pdf, more explanation for the optimization and additional.


## Run parallel code and get binary .dat data result
1. First you need you create two folder in the root path called `result_neu` and `result_per`. These folders will be used to save the parallel data.
2. Go to `parallel_version.cpp`, see the start of main function. There is one line `int bc = 0;`. Set the boundary condition that you want to use, 0 for neumann and 1 for periodic.
3. If you run the code from terminal (set the process number 8 to what you want):
   ```
    mpicxx -o parallel parallel_version.cpp
    mpirun -np 8 ./parallel
   ```
   if you want to submit the code to hpc, first set the parameter in the `evaluate.pbs` to what you want. The default process number is 4. Then:
   ```
    qsub evaluate.pbs
   ```
4. Result Data. If you run with neumann model, the result would be in the folder `result_neu` that you just created. If you run with periodic model, the result would be in the folder `result_per` that you just created. Note these .dat are in binary format. We cannot read it directly.


## Visualization of the parallel data
Parallel data that we just get is not able to be read directly. So we are going to visualize them:
1. Go to `visualization.py`, modify the input path `folder_path` and the output path `output_name1` and `output_name1`. The default input path is folder `result_neu`. The default output paths are `animation/neumann_C1.gif` and `animation/neumann_C2.gif` (__The default output path will override the animation submission folder, please be carefully. You may simply change to different file names for the output path__). 
2. Just run the python code `visualization.py` in any way you like, for example:
   ```
    python3 visualization.py
   ```
   You may need to install the package imageio first:
   ```
    pip install imageio
   ```
3. You will see the visualization gif result of the parallel data in the path that you set. The default path is in the folder animation.


## Extensions and Additions
There are two things that are implemented in this part.
### 1. Fault Tolerance: Checkpoint Mechanism
In this coursework, it is a simple initialize condition and short `max_time`. If we need to run it with more complicated initial condition and far more longer `max_time`. It is very useful to have checkpoint and save the current process regularly. In the `Additional_version.cpp`, the checkpint mechanism is added and you can set the checkpint interval at the beginning of the main function. It will save the current status every time reach the checkpoint. If the execution of the program is interrupted. It can reload the latest checkpoint and recovery from that point when you run the program again. (_Note: Do not forget to create two folders `result_neu` and `result_per` in this directory before execute the code, or you will not able to collect any result data_)

### 2. Testing
One more extension in my work is to write the code to testing the result from parallel code and check whether it is absolute the same with the result data from the serial code.
To run the test with default parameter setting:
- You should have two folders in the root path. One is called `serial_per`, contains the output data from the default serial code, these data are in text format. The other one is called `result_per`, contains the output data from the parallel code, these data are in binary format (_Note those data from both serial and parallel should be generated with the same boundary condition, for example, periodic_).
- If you use different data, modify the setting in the `Test.cpp`. Change `imax` and `jmax` to adapt your data, the default size is 301 * 301. And change the number of iteration condition of the for loop to test with all the data.
- Run the `Test.cpp` with any compiler you like. For example:
  ```
    g++ -o test Test.cpp
    ./test
  ```