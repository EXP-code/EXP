
#include <utility>
#include <deque>
#include <iterator>
#include <algorithm>
#include <vector>
#include <string>

int wordSplit(string& x, vector<string>& words)
{
  words.clear();
  // erase(words.begin(), words.end());

  int i=0;
  while (i<x.size()) {
				// White space
    if (x[i] == ' ' || x[i] == '\t' || x[i] == '\n') {
      i++;
      continue;
    }
				// New word
    string word;
    while (i<x.size() && x[i] != ' ' && x[i] != '\t' && x[i] != '\n')
      word += x[i++];
    word += '\0';
    words.push_back(word);
  }

  return words.size();
}
