#ifndef __JPEG_EXCEPTION__
#define __JPEG_EXCEPTION__

#include <exception>
#include <string>


class FileNotFoundException: public exception {
    FileNotFoundException(std:string message): exception(message) {}
}

#endif