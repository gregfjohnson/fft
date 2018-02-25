-----------------------------------------------------------------------------
-- Copyright (c) Greg Johnson, Gnu Public Licence v. 2.0.
-----------------------------------------------------------------------------
fft = require 'fft'
require 'cmath'
--[[
    Quick-check tests of fft.lua.

    Usage:  lua test_fft.lua

    all tests should be true.
--]]

local function pvec(v) for i = 1, #v do print(v[i]) end end

local function cmp(v1, v2)
    local function eq(x, y) return cmath.abs(x-y) < 1e-15 end

    if #v1 ~= #v2 then return false end

    for i = 1, #v1 do
        if not eq(v1[i], v2[i]) then return false end
    end

    return true
end

print('vector of length 1:  ',         cmp({1},               fft{1}))
print('complex vector of length 1:  ', cmp({1+i},             fft{1+i}))
print('{1,1}:  ',                      cmp({math.sqrt(2), 0}, fft{1,1}))
print('{1,1,1,1}:  ',                  cmp({2,0,0,0},         fft{1,1,1,1}))
print('{1,i,-1,-i}:  ',                cmp({0,2,0,0},         fft{1,i,-1,-i}))

local expected = {
    12.000000000000 + i * 0.000000000000,
    -1.774316085207 + i * 2.258290640619,
     1.289861687269 + i * 0.395365152057,
    -3.000000000000 + i * 1.732050807569,
     1.984454397937 - i * 0.996900084778,
     1.984454397937 + i * 0.996900084778,
    -3.000000000000 - i * 1.732050807569,
     1.289861687269 - i * 0.395365152057,
    -1.774316085207 - i * 2.258290640619,
}
print('{3,1,4,1,5,9,2,6,5}:  ', cmp(expected, fft{3,1,4,1,5,9,2,6,5}))

print('round trip vector of length 1:  ',         cmp({1},                 fft(fft{1}, false)))
print('round trip complex vector of length 1:  ', cmp({1+i},               fft(fft{1+i}, false)))
print('round trip {1,1}:  ',                      cmp({1,1},               fft(fft{1,1}, false)))
print('round trip {1,1,1,1}:  ',                  cmp({1,1,1,1},           fft(fft{1,1,1,1}, false)))
print('round trip {1,i,-1,-i}:  ',                cmp({1,i,-1,-i},         fft(fft{1,i,-1,-i}, false)))
print('round trip {3,1,4,1,5,9,2,6,5}:  ',        cmp({3,1,4,1,5,9,2,6,5}, fft(fft{3,1,4,1,5,9,2,6,5}, false)))
