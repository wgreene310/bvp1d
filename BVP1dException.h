// Copyright (C) 2016 William H. Greene
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

#include <stdexcept>
#include <string>

class BVP1dException : public std::exception{
public:
  BVP1dException(const char *id, const char *msg) :
    id(id), msg(msg) {}
  const char *getId() const { return id.data();  }
#if __cplusplus <= 199711L
  virtual const char *what() const { return msg.data(); }
#else
  virtual const char *what() const noexcept{ return msg.data(); }
#endif
private:
  std::string id, msg;
};