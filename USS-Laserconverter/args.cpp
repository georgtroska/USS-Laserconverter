/*

Copyright (c) 2014 Jarryd Beck

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Visit at https://github.com/jarro2783/cxxopts


*/
#include "stdafx.h"
#include <iostream>
#include <string>

#include "cxxopts.h"
#include "args.h"


void ArgParser::parse(int argc, const char* argv[]) {
    try
    {
        cxxopts::Options options(argv[0], " - example command line options");
        options
            .positional_help("<outputfile> <inputfile> [<inputfile2>,..]");
        //.show_positional_help();

        bool apple = false;

        options
            .set_width(70)
            .set_tab_expansion()
            .allow_unrecognised_options()
            .add_options("Group")
            /*
            ("a,apple", "an apple", cxxopts::value<bool>(apple))
            ("b,bob", "Bob")
            ("char", "A character", cxxopts::value<char>())
            ("t,true", "True", cxxopts::value<bool>()->default_value("true"))
            ("f, file", "File", cxxopts::value<std::vector<std::string>>(), "FILE")
            ("i,input", "Input", cxxopts::value<std::vector<std::string>>())
            ("o,output", "Output file", cxxopts::value<std::string>()
                ->default_value("out.root"))
            ("x", "A short-only option", cxxopts::value<std::string>())
            ("positional",
                "Positional arguments: these are the arguments that are entered "
                "without an option", cxxopts::value<std::vector<std::string>>())
            ("long-description",
                "thisisareallylongwordthattakesupthewholelineandcannotbebrokenataspace")
            ("help", "Print help")
            ("tab-expansion", "Tab\texpansion")
            ("int", "An integer", cxxopts::value<int>(), "N")
            ("float", "A floating point number", cxxopts::value<float>())
            ("vector", "A list of doubles", cxxopts::value<std::vector<double>>())
            ("option_that_is_too_long_for_the_help", "A very long option")
#ifdef CXXOPTS_USE_UNICODE
            ("unicode", u8"A help option with non-ascii: à. Here the size of the"
                " string should be correct")
#endif
            ;



        options.add_options("Group")
            ("c,compile", "compile")
            ("d,drop", "drop", cxxopts::value<std::vector<std::string>>());
*/
        ("help", "Print help")
        ("d,delim", "Delimitor between columns", cxxopts::value<std::string>(_delim1)->default_value("\t"))
        ("v,vecdelim", "Delimitor to sepatate a vector within a column", cxxopts::value<std::string>(_delim2)->default_value(","))
        ("i,input", "Input", cxxopts::value<std::vector<std::string>>(_inputFilenameList))
        ("o,output", "Output file", cxxopts::value<std::string>(_outFilename)
            ->default_value("out.root"))
            ("t,tree", "Name of the tree", cxxopts::value<std::string>(_treename)->default_value("tree"))
            ;
        options.parse_positional({ "input" });

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help({ "","Group" }) << std::endl;

            exit(0);
        }

        if (_inputFilenameList.size() == 0) {
            std::cout << "At least one inputfile needs to be provided!" << std::endl;
            exit(1);
        }

        unsigned int assumeOutfileCount = 0;
        for (unsigned int i = 0; i < _inputFilenameList.size(); i++) {
            std::size_t findstr = _inputFilenameList.at(i).find(".root");
            std::size_t len = _inputFilenameList.at(i).length();
            //Checking if .root is in that string and if it is at the end of the string
            if (findstr != std::string::npos && len - findstr - 1 == 4) {
                if (assumeOutfileCount == 0) {
                    std::cout << "File " << _inputFilenameList.at(i) << " is ending on '.root'. Assuming this as outputfile" << std::endl;
                    _outFilename = _inputFilenameList.at(i);
                    _inputFilenameList.erase(_inputFilenameList.begin() + i);
                    assumeOutfileCount++;
                }
                else {
                    std::cout << "Only one outputfile is allowed!" << std::endl;
                    exit(1);
                }

            }
        }
    }
       
    catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}
