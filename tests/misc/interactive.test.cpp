// Runs a solution built on the snippet in a child process and
// plays the judge through two pipes, as a real interactor does.
#include "../testlib.hpp"
// The Tweak's list printer; say() finds it only if it comes first.
ostream& operator<<(ostream& os, const vector<ll>& v) {
  rep(i, 0, sz(v)) os << (i ? " " : "") << v[i];
  return os;
}
#include "../../content/misc/interactive.hpp"
#include <poll.h>
#include <sys/wait.h>
#include <unistd.h>

// Protocol: T, then per test n. Hidden x in [1, n]; "? y" replies
// 1 iff x <= y; "? y z" replies "y + z" and "y * z" on one line;
// "! x d" is right if d is the first number minus the second (so
// the order of the tokens matters), and the judge confirms with 1.
void solution() {
  // As in the contest template.
  cin.tie(0)->sync_with_stdio(0);
  cin.exceptions(cin.failbit);
  int t = (int)reply();
  while (t--) {
    ll n = reply(), lo = 1, hi = n;
    while (lo < hi) {
      ll mid = lo + (hi - lo) / 2;
      if (ask(mid)) hi = mid;
      else lo = mid + 1;
    }
    ll s = ask(lo, 7), p = reply();
    if (t % 2) answer(lo, s - p);
    else answer(vector<ll>{lo, s - p});
    reply();
  }
}

const int WAIT_MS = 5000; // a missing flush shows up as a timeout

struct Judge {
  int to, from;
  pid_t pid;
  string buf;
  Judge() {
    int a[2], b[2];
    CHECK(pipe(a) == 0 && pipe(b) == 0);
    fflush(stdout), fflush(stderr);
    pid = fork();
    CHECK(pid >= 0);
    if (pid == 0) {
      dup2(a[0], 0), dup2(b[1], 1);
      for (int fd : {a[0], a[1], b[0], b[1]}) close(fd);
      solution();
      exit(0);
    }
    close(a[0]), close(b[1]);
    to = a[1], from = b[0];
  }
  void send(const string& s) {
    string line = s + "\n";
    CHECK(write(to, line.data(), sz(line)) == sz(line));
  }
  // The next line the child wrote, or nullopt on timeout or EOF.
  optional<string> line() {
    while (true) {
      size_t nl = buf.find('\n');
      if (nl != string::npos) {
        string s = buf.substr(0, nl);
        buf.erase(0, nl + 1);
        return s;
      }
      pollfd p{from, POLLIN, 0};
      if (poll(&p, 1, WAIT_MS) <= 0) return nullopt;
      char tmp[4096];
      ssize_t got = read(from, tmp, sizeof tmp);
      if (got <= 0) return nullopt;
      buf.append(tmp, got);
    }
  }
  string fail(const string& msg) {
    kill(pid, SIGKILL);
    finish();
    return msg;
  }
  // Exit code of the child; -1 if it did not exit normally.
  int finish() {
    if (to >= 0) close(to), to = -1;
    int status;
    waitpid(pid, &status, 0);
    close(from);
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
  }
};

// Plays one session; returns an error message, empty if all went
// well. With reject = k >= 0, the judge answers the k-th line the
// child sends with -1 (or, with hang_up, closes the pipe), and the
// child must exit with status 0 without writing more.
string session(const vector<ll>& ns, int reject, bool hang_up, ll& queries) {
  Judge j;
  int lines = 0;
  j.send(to_string(sz(ns)));
  for (ll n : ns) {
    ll x = rnd(1, n), used = 0;
    j.send(to_string(n));
    while (true) {
      optional<string> got = j.line();
      if (!got) return j.fail("no flushed line from the child");
      if (lines++ == reject) {
        if (!hang_up) j.send("-1");
        close(j.to), j.to = -1;
        optional<string> more = j.line();
        if (more) return j.fail("child wrote after -1: " + *more);
        return j.finish() == 0 ? "" : "bad exit after -1";
      }
      istringstream in(*got);
      string kind, rest;
      ll y, z = -1;
      in >> kind >> y;
      bool two = bool(in >> z);
      ostringstream norm;
      norm << kind << ' ' << y;
      if (two) norm << ' ' << z;
      if (norm.str() != *got || (in >> rest))
        return j.fail("malformed line: " + *got);
      if (kind == "?" && !two) {
        if (y < 1 || y > n) return j.fail("query out of range");
        if (++used > __lg(2 * n - 1)) return j.fail("too many");
        j.send(x <= y ? "1" : "0");
      } else if (kind == "?") {
        j.send(to_string(y + z) + " " + to_string(y * z));
      } else if (kind == "!" && two) {
        if (y != x || z != (x + 7) - x * 7)
          return j.fail("wrong answer " + *got);
        j.send("1");
        break;
      } else return j.fail("unknown line " + *got);
    }
    queries += used + 2;
  }
  if (optional<string> more = j.line()) return j.fail("extra output: " + *more);
  return j.finish() == 0 ? "" : "bad exit status";
}

int main() {
  signal(SIGPIPE, SIG_IGN);
  ll queries = 0;
  rep(it, 0, iters(150)) {
    vector<ll> ns(rnd(1, 4));
    for (ll& n : ns) n = rnd(0, 1) ? rnd(1, 20) : rnd(1, 1e18);
    int reject = rnd(0, 2) ? -1 : (int)rnd(0, 6);
    CHECK_EQ(session(ns, reject, rnd(0, 1), queries), string());
  }
  // Speed: 3000 tests of 32 queries each, about 1e5 round trips.
  double start = elapsed();
  ll before = queries;
  CHECK_EQ(session(vector<ll>(3000, (ll)1e9), -1, 0, queries), string());
  fprintf(stderr, "  %lld round trips: %.2fs\n", queries - before,
          elapsed() - start);
  puts("OK misc/interactive");
}
