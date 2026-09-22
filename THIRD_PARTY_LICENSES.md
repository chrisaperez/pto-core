# Third-party licenses

`pto-core` itself is MIT-licensed (see [`LICENSE`](LICENSE)). It vendors three
third-party libraries, each under its own license, listed here with the
copyright holder, version, and location in this tree. None of these licenses
require anything beyond what this file and the vendored copies already do:
MIT and Apache-2.0 both permit redistribution as long as the notice travels
with the code, which is what this file is for.

| Library | Location | License | Copyright |
|---|---|---|---|
| [cpp-httplib](https://github.com/yhirose/cpp-httplib) | `modules/cuttag_profiler/third_party/httplib/` | MIT | (c) 2024 Yuji Hirose |
| [nlohmann/json](https://github.com/nlohmann/json) | `modules/cuttag_profiler/third_party/nlohmann/` | MIT | (c) 2013-2023 Niels Lohmann |
| [hnswlib](https://github.com/nmslib/hnswlib) | `modules/scrna_matrix/third_party/hnswlib/` | Apache-2.0 | hnswlib contributors |

`hnswlib` ships its own license file in-tree at
[`modules/scrna_matrix/third_party/hnswlib/LICENSE`](modules/scrna_matrix/third_party/hnswlib/LICENSE);
it is not reproduced below. It is vendored at v0.9.0 with four local
concurrency/memory-safety patches applied on top — see that directory's
`VERSION.txt` for exactly what changed and why. `cpp-httplib` and
`nlohmann/json` are vendored as single unmodified headers and carry their
license only as the SPDX/copyright comment at the top of the file, so the
full MIT text for each is reproduced here.

---

## cpp-httplib — MIT License

```
Copyright (c) 2024 Yuji Hirose. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## nlohmann/json — MIT License

```
SPDX-FileCopyrightText: 2013-2023 Niels Lohmann <https://nlohmann.me>
SPDX-License-Identifier: MIT

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## hnswlib — Apache License 2.0

See [`modules/scrna_matrix/third_party/hnswlib/LICENSE`](modules/scrna_matrix/third_party/hnswlib/LICENSE)
for the full text. Summary only: Apache-2.0 permits use, modification and
redistribution (including in a binary/wheel), requires the license and any
`NOTICE` file to travel with copies, and requires stating that modified files
were changed — which `modules/scrna_matrix/third_party/hnswlib/VERSION.txt`
already does for this vendored copy's four local patches.
