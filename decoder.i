%module decoder
%{ 
    #include "JPEGDecoder.hpp"
%}

%rename(_print) operator <<;