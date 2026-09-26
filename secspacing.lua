-- pandoc-crossref (cref:false) separates the section symbol from the reference
-- number with a non-breaking space, printing "§ 4.2"; this normalizes it to "§4.2"
-- by rewriting the LaTeX raw output, after -F pandoc-crossref has run.

local SECTION = "\194\167" -- U+00A7 SECTION SIGN
local NBSP    = "\194\160" -- U+00A0 NO-BREAK SPACE

local function tighten(s)
  if type(s) ~= "string" then return s end
  s = s:gsub(SECTION .. "~", SECTION)
  s = s:gsub(SECTION .. NBSP, SECTION)
  s = s:gsub(SECTION .. " ", SECTION)
  return s
end

local function handleRaw(el)
  if el.format == "latex" or el.format == "tex" then
    el.text = tighten(el.text)
  end
  return el
end

local function handleStr(el)
  el.text = tighten(el.text)
  return el
end

return {
  { Str = handleStr, RawInline = handleRaw, RawBlock = handleRaw },
}
