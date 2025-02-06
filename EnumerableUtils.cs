#nullable disable

namespace LinearAlgebra
{
    public static class EnumerableUtils
    {
        public static T Map<T>(this IEnumerable<T> source, Func<T, T, T> map)
        {
            T res = default;
            foreach (var item in source)
                res = map(res, item);
            return res;
        }

        public static string RemoveSpaces(this string s)
        {
            var charArray = s.ToCharArray();
            charArray = charArray.Where(x => !char.IsWhiteSpace(x)).ToArray();

            return new(charArray);
        }
    }
}
