// Example for C++ Interface to Gnuplot

// requirements:
// * gnuplot has to be installed (http://www.gnuplot.info/download.html)
// * for Windows: set Path-Variable for Gnuplot path (e.g. C:/program files/gnuplot/bin)
//             or set Gnuplot path with: Gnuplot::set_GNUPlotPath(const std::string &path);


#include <iostream>
#include "gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communikation with Gnuplot

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
 #include <conio.h>   //for getch(), needed in wait_for_key()
 #include <windows.h> //for Sleep()
 void sleep(int i) { Sleep(i*1000); }
#endif


#define SLEEP_LGTH 2  // sleep time in seconds
#define NPOINTS    50 // length of array

void wait_for_key(); // Programm halts until keypress

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
    // if path-variable for gnuplot is not set, do it with:
    // Gnuplot::set_GNUPlotPath("C:/program files/gnuplot/bin/");

    // set a special standard terminal for showonscreen (normally not needed),
    //   e.g. Mac users who want to use x11 instead of aqua terminal:
    // Gnuplot::set_terminal_std("x11");

    cout << "*** example of gnuplot control through C++ ***" << endl << endl;

    //
    // Using the GnuplotException class
    //
    try
    {
        Gnuplot g1("lines");

   


        //
        // User defined 1d, 2d and 3d point sets
        //
        std::vector<double> x, y, y2, dy, z;

        for (int i = 0; i < NPOINTS; i++)  // fill double arrays x, y, z
        {
            x.push_back((double)i);             // x[i] = i
            y.push_back((double)i * (double)i); // y[i] = i^2
            z.push_back( x[i]*y[i] );           // z[i] = x[i]*y[i] = i^3
            dy.push_back((double)i * (double)i / (double) 10); // dy[i] = i^2 / 10
        }
//         y2.push_back(0.00); y2.push_back(0.78); y2.push_back(0.97); y2.push_back(0.43);
//         y2.push_back(-0.44); y2.push_back(-0.98); y2.push_back(-0.77); y2.push_back(0.02);


        g1.reset_all();
//         cout << endl << endl << "*** user-defined lists of doubles" << endl;
//         g1.set_style("impulses").plot_x(y,"user-defined doubles");

        g1.reset_plot();
        cout << endl << endl << "*** user-defined lists of points (x,y)" << endl;
        g1.set_grid();
        g1.set_style("points").plot_xy(x,z,"user-defined points 2d");

  

  


        //
        // Multiple output screens
        //
        cout << endl << endl;
        cout << "*** multiple output windows" << endl;

    

      

        

        wait_for_key();

    }
    catch (GnuplotException ge)
    {
        cout << ge.what() << endl;
    }


    cout << endl << "*** end of gnuplot example" << endl;

    return 0;
}



void wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}
