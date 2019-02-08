////////////////////////////////////////////
//
// A C++ interface to gnuplot. 
//
// This is a direct translation from the C interface
// written by N. Devillard (which is available from
// http://ndevilla.free.fr/gnuplot/).
//
// As in the C interface this uses pipes and so wont
// run on a system that does'nt have POSIX pipe 
// support
//
// Rajarshi Guha
// <rajarshi@presidency.com>
//
// 07/03/03
//
// /////////////////////////////////////////

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

#include <stdarg.h>
#include <unistd.h>

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include <stdexcept>

#define GP_MAX_TMP_FILES    64
#define GP_TMP_NAME_SIZE    512
#define GP_CMD_SIZE         1024
#define GP_TITLE_SIZE       80

#define PATH_MAXNAMESZ       4096

using namespace std;

template <typename Container>
void
	stringtok (Container &container, string const &in,
const char * const delimiters = " \t\n")
{
	const string::size_type len = in.length();
	string::size_type i = 0;

	while ( i < len )
	{
		// eat leading whitespace
		i = in.find_first_not_of (delimiters, i);
		if (i == string::npos)
			return;   // nothing left but white space

		// find the end of the token
		string::size_type j = in.find_first_of (delimiters, i);

		// push token
		if (j == string::npos) {
			container.push_back (in.substr(i));
			return;
		} else
			container.push_back (in.substr(i, j-i));

		// set up for next loop
		i = j + 1;
	}
}

using namespace std;

class GnuplotException : public runtime_error
{
public:
	GnuplotException(const string &msg) : runtime_error(msg){}
};

class Gnuplot
{
private:
	
	FILE            *gnucmd;
	string           pstyle;
	vector<string>   to_delete;
	int              nplots;

	bool get_program_path(const string pname)
	{
		list<string> ls;
		char *path;

		path = getenv("PATH");
		if (!path)
		{
			cerr << "Path is not set" << endl;
			return false;
		}
		else
		{
			stringtok(ls,path,":");
			for (list<string>::const_iterator i = ls.begin();
			i != ls.end(); ++i)
			{
				string tmp = (*i) + "/" + pname;
				if (access(tmp.c_str(),X_OK) == 0)
					return true;
			}
		}
		return false;
	}
		
public:
		
	Gnuplot() 
	{
		if (getenv("DISPLAY") == NULL)
			throw GnuplotException("cannot find DISPLAY variable");
		if (!this->get_program_path("gnuplot"))
			throw GnuplotException("Can't find gnuplot in your PATH");
    
		this->gnucmd = popen("gnuplot","w");
		if (!this->gnucmd)
			throw GnuplotException("Could'nt open connection to gnuplot");

		this->set_style("points");
		this->nplots = 0;
	}

	// set a style during construction
	Gnuplot(const string & style) 
	{
		if (getenv("DISPLAY") == NULL)
			throw GnuplotException("cannot find DISPLAY variable");
		if (!this->get_program_path("gnuplot"))
			throw GnuplotException("Can't find gnuplot in your PATH");

		this->gnucmd = popen("gnuplot","w");
		if (!this->gnucmd)
			throw GnuplotException("Could'nt open connection to gnuplot");
		this->set_style(style);
		this->nplots = 0;
	}
        
	// The equivilant of gnuplot_plot_once, the two forms
	// allow you to plot either (x,y) pairs or just a single
	// vector at one go
	Gnuplot(const string &title, const string &style, const string &labelx,  const string &labely,
	vector<double> x, vector<double> y)
	{
		if (!this->get_program_path("gnuplot"))
			throw GnuplotException("Can't find gnuplot in your PATH");

		if (getenv("DISPLAY") == NULL)
			throw GnuplotException("cannot find DISPLAY variable");

		this->gnucmd = popen("gnuplot","w");
		if (!this->gnucmd)
			throw GnuplotException("Could'nt open connection to gnuplot");
		this->nplots = 0;


		if (x.size() == 0 || y.size() == 0)
			throw GnuplotException("vectors too small");
		if (!this->gnucmd)
			throw GnuplotException("Could'nt open connection to gnuplot");

		if (style == "")
			this->set_style("lines");
		else
			this->set_style(style);

		if (labelx == "")
			this->set_xlabel("X");
		else
			this->set_xlabel(labelx);
		if (labely == "")
			this->set_ylabel("Y");
		else
			this->set_ylabel(labely);
    
		this->plot_xy(x,y,title);

		cout << "Press enter to continue" << endl;
		while (getchar() != '\n'){}
	}
        
	Gnuplot(const string &title, const string &style, const string &labelx,  const string &labely, vector<double> x)
	{
		if (!this->get_program_path("gnuplot"))
			throw GnuplotException("Can't find gnuplot in your PATH");

		if (getenv("DISPLAY") == NULL)
			throw GnuplotException("cannot find DISPLAY variable");

		this->gnucmd = popen("gnuplot","w");
		if (!this->gnucmd)
			throw GnuplotException("Could'nt open connection to gnuplot");
		this->nplots = 0;


		if (x.size() == 0)
			throw GnuplotException("vector too small");
		if (!this->gnucmd)
			throw GnuplotException("Could'nt open connection to gnuplot");

		if (style == "")
			this->set_style("lines");
		else
			this->set_style(style);

		if (labelx == "")
			this->set_xlabel("X");
		else
			this->set_xlabel(labelx);
		if (labely == "")
			this->set_ylabel("Y");
		else
			this->set_ylabel(labely);
    
		this->plot_x(x,title);

		cout << "Press enter to continue" << endl;
		while (getchar() != '\n'){}
	}
        
	~Gnuplot() 
	{
		if ((this->to_delete).size() > 0)
		{
			for (int i = 0; i < this->to_delete.size(); i++)
				remove(this->to_delete[i].c_str());
		}
		if (pclose(this->gnucmd) == -1)
			cerr << "Problem closing communication to gnuplot" << endl;
		return;
	}

	// send a command to gnuplot
	void cmd(const char * cmdstr, ...)
	{
		va_list ap;
		char local_cmd[GP_CMD_SIZE];

		va_start(ap, cmdstr);
		vsprintf(local_cmd, cmdstr, ap);
		va_end(ap);
		strcat(local_cmd,"\n");
		fputs(local_cmd,this->gnucmd);
		fflush(this->gnucmd);
		return;
	}

	// set line style
	void set_style(const string &stylestr)
	{
		if (stylestr != "lines" &&
			stylestr != "points" &&
			stylestr != "linespoints" &&
			stylestr != "impulses" &&
			stylestr != "dots" &&
			stylestr != "steps" &&
			stylestr != "errorbars" &&
			stylestr != "boxes" &&
			stylestr != "boxerrorbars")
			this->pstyle = string("points");
		else
			this->pstyle = stylestr;
	}

	// set y and x axis labels
	void set_ylabel(const string &label)
	{
		ostringstream cmdstr;

		cmdstr << "set xlabel \"" << label << "\"";
		auto css = cmdstr.str();
		const char* cstr2 = css.c_str(); 
		this->cmd(cstr2);

		return;
	}

	void set_xlabel(const string &label)
	{
		ostringstream cmdstr;

		cmdstr << "set xlabel \"" << label << "\"";
		auto css = cmdstr.str();
		const char* cstr2 = css.c_str(); 
		this->cmd(cstr2);

		return;
	}

	// plot a single vector
	void plot_x(vector<double> d, const string &title)
	{
		ofstream tmp;
		ostringstream cmdstr;
		char name[] = "/tmp/gnuplotiXXXXXX";

		if (this->to_delete.size() == GP_MAX_TMP_FILES - 1)
		{
			cerr << "Maximum number of temporary files reached (" << GP_MAX_TMP_FILES << "): cannot open more files" << endl;
			return;
		}

		//
		//open temporary files for output
		if (mkstemp(name) == -1)
		{
			cerr << "Cannot create temporary file: exiting plot" << endl;
			return;
		}
		tmp.open(name);
		if (tmp.bad())
		{
			cerr << "Cannot create temorary file: exiting plot" << endl;
			return;
		}

		//
		// Save the temporary filename
		// 
		this->to_delete.push_back(name);

		//
		// write the data to file
		//
		for (int i = 0; i < d.size(); i++)
			tmp << d[i] << endl;
		tmp.flush();    
		tmp.close();

		//
		// command to be sent to gnuplot
		//
		if (this->nplots > 0)
			cmdstr << "replot ";
		else cmdstr << "plot ";
		if (title == "")
			cmdstr << "\"" << name << "\" with " << this->pstyle;
		else
			cmdstr << "\"" << name << "\" title \"" << title << "\" with " << this->pstyle;

		//
		// Do the actual plot
		//
		auto css = cmdstr.str();
		const char* cstr2 = css.c_str(); 
		this->cmd(cstr2);
		this->nplots++;

		return;
	}

	// plot x,y pairs
	void plot_xy(vector<double> x, vector<double> y, const string &title)
	{
		ofstream tmp;
		ostringstream cmdstr;
		char name[] = "/tmp/gnuplotiXXXXXX";
    
		// should raise an exception
		if (x.size() != x.size())
			return;

		if ((this->to_delete).size() == GP_MAX_TMP_FILES - 1)
		{
			cerr << "Maximum number of temporary files reached (" << GP_MAX_TMP_FILES << "): cannot open more files" << endl;
			return;
		}

		//
		//open temporary files for output
		//
		if (mkstemp(name) == -1)
		{
			cerr << "Cannot create temporary file: exiting plot" << endl;
			return;
		}
		tmp.open(name);
		if (tmp.bad())
		{
			cerr << "Cannot create temorary file: exiting plot" << endl;
			return;
		}

		//
		// Save the temporary filename
		// 
		this->to_delete.push_back(name);

		//
		// write the data to file
		//
		for (int i = 0; i < x.size(); i++)
			tmp << x[i] << " " << y[i] << endl;
		tmp.flush();    
		tmp.close();

		//
		// command to be sent to gnuplot
		//
		if (this->nplots > 0)
			cmdstr << "replot ";
		else cmdstr << "plot ";
		if (title == "")
			cmdstr << "\"" << name << "\" with " << this->pstyle;
		else
			cmdstr << "\"" << name << "\" title \"" << title << "\" with " << this->pstyle;

		//
		// Do the actual plot
		//
		auto css = cmdstr.str();
		const char* cstr2 = css.c_str(); 
		this->cmd(cstr2);
		this->nplots++;

		return;
	}

	// plot an equation of the form: y = ax + b
	// You supply a and b
	void plot_slope(double a, double b, const string &title)
	{
		ostringstream stitle;
		ostringstream cmdstr;

		if (title == "")
			stitle << "no title";
		else
			stitle << title;

		if (this->nplots > 0)
			cmdstr << "replot " << a << " * x + " << b << " title \"" << stitle.str() << "\" with " << pstyle;
		else
			cmdstr << "plot " << a << " * x + " << b << " title \"" << stitle.str() << "\" with " << pstyle;
		auto css = cmdstr.str();
		const char* cstr2 = css.c_str(); 
		this->cmd(cstr2);
		this->nplots++;
		return;
	}

	// plot an equation supplied as a string
	void plot_equation(const string &equation, const string &title)
	{
		string titlestr, plotstr;
		ostringstream cmdstr;

		if (title == "")
			titlestr = "no title";
		else
			titlestr = title;

		if (this->nplots > 0)
			plotstr = "replot";
		else
			plotstr = "plot";

		cmdstr << plotstr << " " << equation << " " << "title \"" << titlestr << "\" with " << this->pstyle;
		auto css = cmdstr.str();
		const char* cstr2 = css.c_str(); 
		this->cmd(cstr2);
		this->nplots++;

		return;
	}
	
	// if multiple plots are present it will clear the plot area
	void reset_plot(void)
	{       
		if (this->to_delete.size() > 0)
		{
			for (int i = 0; i < this->to_delete.size(); i++)
				remove(this->to_delete[i].c_str());
		}
		this->nplots = 0;
		return;
	}
        
};

#endif
