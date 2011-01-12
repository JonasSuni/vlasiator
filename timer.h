#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <vector>

#ifdef PROFILE
/** A simple class for timing the wall-clock time spent in 
 * functions. Using Timer class one can create, start and stop 
 * timers. The timers are based on std::clock().
 */
class Timer {
 public:
   Timer();
   ~Timer();
   
   static unsigned int create(const std::string& name);
   static double getValue(const unsigned int& timerID);
   static void print();
   static void start(const unsigned int& timerID);
   static void stop(const unsigned int& timerID);
   
 private:
   /** Struct TimerData contains the variables necessary for 
    * measuring the total time spent between calls to Timer::start 
    * and Timer::stop member functions.
    */
   struct TimerData {
      std::string name; /**< The name of the timer.*/
      std::clock_t startClock; /**< The value of clock() when the timer was started via call to Timer::start.*/
      double timeInSeconds; /**< The total time in seconds for the timer.*/
   };
   
   static std::vector<TimerData> timers; /**< Container for all timers.*/
};

#else 

class Timer {
 public:
   Timer() { }
   ~Timer() { }
   
   static unsigned int create(const std::string& name) {return 0;}
   static double getValue(const unsigned int& timerID) {return 0.0;}
   static void print() { }
   static void start(const unsigned int& timerID) { }
   static void stop(const unsigned int& timerID) { }
 private:
   
};

#endif // #ifdef PROFILE

#endif
