#include <iostream>
#include <cstdlib>
#include <string>
#include <list>
#include <algorithm>
 
class ConsoleMenu {
public:
    // Generic prompt for input
    ConsoleMenu(): msg("Select an option: ") { }
    // Custom prompt, just in case
    ConsoleMenu(std::string omsg): msg(omsg) { }
 
    // The order of 'values' doesn't matter, 'options' does.
    void add(int val, std::string opt)
    {
        values.push_back(val);
        options.push_back(opt);
    }
 
    // For switch statements and the like
    int  opt() { return option; }
    // Interactive input
    bool selection();
 
    // Return ostream& for linked output (e.g. obj.display()<<endl;)
    std::ostream& display();
	
private:
    int                    option;  // Interactive input option
    std::string            msg;     // Prompt message for input
    std::list<std::string> options; // List of printable options
    std::list<int>         values;  // List of valid option values
};
 
std::ostream& ConsoleMenu::display()
{
    std::list<std::string>::const_iterator it = options.begin();
 
    while (it != options.end())
        std::cout << *it++ << '\n';
 
    // Be sure to flush the stream
    return std::cout << msg << std::flush;
}
 
bool ConsoleMenu::selection()
{
    if (std::cin >> option &&
        std::find(values.begin(), values.end(), option) != values.end()) {
        return true; // Good
    }
     
    // cin's cleanup is ugly
    if (!std::cin.good()) {
        std::cin.clear();
         
        int ch;
        while ((ch = std::cin.get()) != '\n' && ch != EOF)
            ;
    }
 
    return false; // Bad
}
 
using namespace std;
 
void a() { cout << "\nOption 1\n" << endl; }
void b() { cout << "\nOption 2\n" << endl; }
void c() { exit(0); }
 
int main()
{
    ConsoleMenu menu;
    void (*action[])() = {a, b, c};
 
    menu.add(1, "1) Option 1");
    menu.add(2, "2) Option 2");
    menu.add(3, "3) Quit");
 
    while (1)
    {
        menu.display();
 
        if (menu.selection())
            action[menu.opt()-1]();
        else
            cerr<<"\nInvalid option\n"<<endl;
    }
}