// consoleUtil.hpp
// a collection of useful interface functions
// for a C++11 console application

#ifndef CONSOLEUTIL_CPP
#define CONSOLEUTIL_CPP

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <functional>

// Clear the screen
void cls() {
  // Windows:
  //system("cls");
    
  // Linux/Unix:
  system("clear");
  //std::cout << "\x1B[2J\x1B[H";
}

// Print a newline/blank line to console
void newl(int n = 1) {
  for (int i = 0; i < n; i++)
    std::cout << std::endl;
}

// Return a string that allows a sentence/spaces
std::string get_s() {
  std::string input;
  getline(std::cin, input);
  return input;
}

// Get a number of type T
template <class T>
T get_n(std::string prompt) {
  T number;
  std::string input;
  while (true) {
    std::cout << prompt;
    input = get_s();
    // converts from string to number safely
    std::stringstream myStream(input);
    if (myStream >> number)
      break;
    else
      std::cout << input << " -> Invalid input" << std::endl;
  }
  return number;
}

// Ask user for a number within a specified range
template<class T>
T get_n_range(std::string prompt, T low, T high) {
  T number;
  std::string input;
  while (true) {
    std::cout << prompt << " [" << low << " - " << high << "]: ";
    input = get_s();
    // converts from string to number safely
    std::stringstream myStream(input);
    if (myStream >> number) {
      if (number >= low && number <= high)
        break;
    }
    std::cout << input << " -> Invalid input" << std::endl;
  }
  return number;
}

// Ask user for a Y/N reply to a question
bool yesNo(std::string question) {
  // function returns true for yes, false for no
  char reply = { 0 };
  std::string input;
  while (true) {
    std::cout << question << " (y/n): ";
    input = get_s();
    if (input.length() == 1)
    {
      reply = toupper(input[0]);
      if (reply == 'Y' || reply == 'N')
        break;
    }
    std::cout << input << " -> Invalid input." << std::endl;
  }
  return (reply == 'Y' ? true : false);
}

// Print a waiting message. Wait for Enter/Return
void waitReturn(std::string message) {
  std::cout << message << " ";
  get_s();
}

// print a title with underline or border
void title(std::string title, char ch, bool border) {
  if (!border) {
    std::cout << title << std::endl;
    std::cout << std::string(title.length(), ch) << std::endl;
  }
  else {
    // with a border
    std::cout << std::string(title.length() + 8, ch) << std::endl;
    std::cout << std::string(2, ch) << "  " << std::string(title.length(), ' ') 
      << "  " << std::string(2, ch) << std::endl;
    std::cout << std::string(2, ch) << "  " << title << "  " << std::string(2, ch) << std::endl;
    std::cout << std::string(2, ch) << "  " << std::string(title.length(), ' ') 
      << "  " << std::string(2, ch) << std::endl;
    std::cout << std::string(title.length() + 8, ch) << std::endl;
  }
}

// show a menu from a list of items
void menuList(const std::vector<std::string>& items) {
  for (unsigned int i = 0; i < items.size(); i++) {
    std::cout << "  " << items[i] << '\n';
  }
  std::cout << std::endl;
}

int menu(std::string mytitle, const std::vector<std::pair<std::string, std::function<void()> > >& items) {
  title(mytitle, '=', false);
  for (size_t i = 0; i < items.size(); i++) {
    std::cout << "  " << i + 1 << " " << items[i].first << '\n';
  }
  int ans = get_n_range<int>(">", 1, items.size());
  if (items[ans - 1].second != nullptr) items[ans - 1].second();
  return ans;
}

#if 1

void hello();

void mainMenu() {
  int choice = -1;
  while (choice != 4) {
    newl(); 
    choice = menu("Main Menu", {
      {"Nothing", nullptr},
      {"Classic function", hello},
      {"Lambda function", [](){ std::cout << "I'm a Lamdba\n"; }},
      {"Quit", nullptr}
    });
  }
}

void hello() {
   newl(); 
   int choice = menu("Hello Menu", {
     {"Nothing", nullptr},
     {"Say Hello", [](){ std::cout << "Hello!!\n"; }},
     {"Back", mainMenu}
   });
}

// test functions/usage:
int main()
{
  mainMenu();
  return 0;
  //std::cout << "You choice is " << choice << "\n";
  
  /*
  bool lee = yesNo("Is it you Man?");
  if (lee)
    std::cout << "Hello Man!\n";
  else
    std::cout << "Hello Somebody.\n";
  
  int age = get_n_range<int>("How old are you?", 0, 100);
  if (age <= 33)
    std::cout << "Young boy!\n";
  else
    std::cout << "You are so old! (" << age << ")\n";
  
  double fave = get_n<double>("Enter a float number: ");
  std::cout << "You entered " << fave;
  newl(2);
  waitReturn("<Hit RETURN to continue>");
  newl();
  return 0;
  */
}
#endif

#endif /* end of include guard: CONSOLEUTIL_CPP */
