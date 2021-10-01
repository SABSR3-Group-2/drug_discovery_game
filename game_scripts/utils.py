import textwrap

def split_text(text, width=75):
    """Given a block of text and a specified number of characters to display per line, splits the text onto different
    lines.
    :param text: text
    :type text: str
    :param width: number of characters per line
    :type width: int
    :return: text split into lines, stored as a list
    :rtype: list
    """
    length = len(text)
    num_lines = int(length / width)  # number of full lines
    remainder = length % width  # number of characters on non-full line
    paragraph = []
    idx = 0
    for i in range(1, num_lines + 1):
        paragraph.append(text[idx:i*width])
        idx += width
    if remainder != 0:
        paragraph.append(text[-remainder:])
    return paragraph

text = 'hello my name is James and this is some test text. I guess there will be some formatting issues with the lines' \
       'running over and I"ll probably need to include some rules for how to split text'

z = textwrap.fill(text, 20)
print(z)

print(len(text))
print(split_text(text, 20))