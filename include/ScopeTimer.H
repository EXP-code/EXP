#pragma once

#include <chrono>
#include <iostream>
#include <string>

class ScopeTimer {
public:
  using Clock = std::chrono::steady_clock;

  // Define units for output
  enum class Unit {
    Seconds,
    Milliseconds,
    Microseconds,
    Nanoseconds
  };

  // Default to Seconds if the user doesn't choose
  explicit ScopeTimer(std::string task_name, Unit output_unit = Unit::Seconds) 
    : name(std::move(task_name)), unit(output_unit), start_time(Clock::now()) {}
  
  // Output on deletion of the object
  ~ScopeTimer()
  {
    const auto end_time = Clock::now();
    const auto delta = end_time - start_time;
    
    // Convert natively using std::chrono based on user choice
    switch (unit) {
    case Unit::Seconds: {
      std::chrono::duration<double> elapsed = delta;
      std::cout << "[" << name << "] completed in " << elapsed.count() << " s\n";
      break;
    }
    case Unit::Milliseconds: {
      std::chrono::duration<double, std::milli> elapsed = delta;
      std::cout << "[" << name << "] completed in " << elapsed.count() << " ms\n";
      break;
    }
    case Unit::Microseconds: {
      std::chrono::duration<double, std::micro> elapsed = delta;
      std::cout << "[" << name << "] completed in " << elapsed.count() << " µs\n";
      break;
    }
    case Unit::Nanoseconds: {
      std::chrono::duration<double, std::nano> elapsed = delta;
      std::cout << "[" << name << "] completed in " << elapsed.count() << " ns\n";
      break;
    }
    default:
      std::cerr << "Unknown time unit for ScopeTimer\n";
      break;
    }
  }

  ScopeTimer(const ScopeTimer&) = delete;
  ScopeTimer& operator=(const ScopeTimer&) = delete;

private:
  std::string name;
  Unit unit;
  std::chrono::time_point<Clock> start_time;
};
