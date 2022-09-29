#include "LoKI-MC/Headers/Message.h"
#include <iostream>
#include <cstdio>
#include <string>

void Message::error(std::string message){
	// write first in the errorLog file
	FILE* fileID = std::fopen("errorLog.txt", "a");
	std::fprintf(fileID, "Program stopped due to the following error:\n%s\n", message.c_str());
	std::fclose(fileID);
	// write in the terminal
	setTerminalColor("red");
	std::printf("Program stopped due to the following error:\n%s\n", message.c_str());
	setTerminalColor("reset");
	std::exit(EXIT_FAILURE);
}
void Message::warning(std::string message){
	// write first in the errorLog file
	FILE* fileID = std::fopen("errorLog.txt", "a");
	std::fprintf(fileID, "Pay attention to the following warning:\n%s\n", message.c_str());
	std::fclose(fileID);
	setTerminalColor("boldYellow");
	// write in the terminal
	std::printf("Pay attention to the following warning:\n%s\n", message.c_str());
	setTerminalColor("reset");
}
void Message::setTerminalColor(std::string color){
	if (color == "red"){
		std::printf("\033[31m");
	}
	else if (color == "boldYellow"){
		std::printf("\033[1;33m");
	}
	else if (color == "reset"){
		std::printf("\033[0m");
	}
	else{
		std::printf("\033[31m");
		std::printf("The color '%s' used in the argument of setTerminalColor is not valid. Please check the file 'TerminalColors.h'.\n", color.c_str());
		std::printf("\033[0m)");
		std::exit(EXIT_FAILURE);
	}
}