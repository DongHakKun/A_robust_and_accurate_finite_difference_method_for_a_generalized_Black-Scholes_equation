#pragma once
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <xcpp/xdisplay.hpp>

inline void displaySVG(const char* path)
{
    std::ifstream file(path);
    if(!file) throw std::runtime_error("Cannot read SVG figure");
    std::string svg((std::istreambuf_iterator<char>(file)),
                    std::istreambuf_iterator<char>());
    nlohmann::json data;
    data["image/svg+xml"] = svg;
    xeus::get_interpreter().display_data(data, nlohmann::json::object(),
                                        nlohmann::json::object());
}
