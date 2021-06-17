import lz4.frame
import base64


def compress_sequence(sequence):
    return base64.b64encode(lz4.frame.compress(sequence.encode('utf-8'))).decode()


def decompress_sequence(s):
    return lz4.frame.decompress(base64.b64decode(s)).decode('utf-8')

