#ifndef __Message__
#define __Message__
#include <string>

namespace Message{
	void error(std::string message);
	void warning(std::string message);
	void setTerminalColor(std::string color);
};

#endif