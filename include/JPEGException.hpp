#ifndef __JPEG_EXCEPTION__
#define __JPEG_EXCEPTION__

#include <exception>
#include <sstream>
#include <string>


class FileNotFoundException: public std::exception {
    public:
        FileNotFoundException(std::string file): file_(file) {}

        virtual const char* what() const throw() {
            std::stringstream stream;
            stream << "File " << this->file_ << " not found.";
            return stream.str().c_str();
        }

    private:
        std::string file_;
};

#endif