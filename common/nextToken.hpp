#ifndef NEXTTOKEN_HPP_29CD58E6
#define NEXTTOKEN_HPP_29CD58E6

#include <string>
#include <iostream>
#include <vector>
#include <cctype>

// Usage:
// std::string token;
// nextToken tokenizer(stream);
// tokenizer.getNext(token);
class nextToken
{
public:
	std::istream & is;
	std::vector<char> ignored;
	char targetCharacter;
	bool reachTarget;
	
	nextToken(std::istream & IS) : is(IS) {
		// Default tokens that are ignored (excepted space characters)
		ignored.push_back('=');
		ignored.push_back('(');
		ignored.push_back(')');
		ignored.push_back(',');
		ignored.push_back(';');
		ignored.push_back('\"');
		ignored.push_back('{');
		ignored.push_back('}');
		ignored.push_back('[');
		ignored.push_back(']');

		targetCharacter = '}';
		reachTarget = false;
	}

	~nextToken() { ignored.clear(); }


	void setTarget(char C) {
		targetCharacter = C;
		reachTarget = false;
	}

	bool isIgnored(char C) {
		if (isspace(C)) return true;
		if (C == targetCharacter) {
			reachTarget = true;
			return false;
		}
		for (size_t i = 0 ; i < ignored.size() ; i++) {
			if (C == ignored[i]) return true;
		}
		return false;
	}


	void getNext(std::string & token) {
		token.clear();
		char C;

		// Ignored characters
		while (is.good()) {
			C = (char)(is.get());
			if (C == 34) {
				// Read up to the next quotation mark
				C = (char)(is.get());
				while (is.good() && C != 34) {
					token += C;
					C = (char)(is.get());
				}
				return;
			}
			if (isIgnored(C) == false) {
				token += C;
				break;
			}
		}

		// Now read the token up to the next ignored character
		while (is.good()) {
			C = (char)(is.get());
			if (isIgnored(C) == true) break;
			else {
				token += C;
			}
		}

		// Seek to the next readable token (with >>)
		while (is.good()) {
			if (isIgnored(C) == false) {
				is.seekg(-1, std::ios::cur);
				return;
			}
			C = (char)(is.get());
		}
	}
};


#endif /* end of include guard: NEXTTOKEN_HPP_29CD58E6 */
