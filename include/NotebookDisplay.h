#pragma once
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <xcpp/xdisplay.hpp>
#include <xeus/xbase64.hpp>

inline void displayPNG(const char* path)
{
    std::ifstream file(path, std::ios::binary);
    if(!file) throw std::runtime_error("Cannot read PNG figure");
    std::string png((std::istreambuf_iterator<char>(file)),
                    std::istreambuf_iterator<char>());
    nlohmann::json data;
    data["image/png"] = xeus::base64encode(png);
    xeus::get_interpreter().display_data(data, nlohmann::json::object(),
                                        nlohmann::json::object());
}
