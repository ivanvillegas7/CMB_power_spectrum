# -*- coding: utf-8 -*-
import requests
import io
import base64

import matplotlib.pyplot as plt

def bot_texter(bot_message, file_data=None, file_name=None, caption=None):

    """
    Send a notification through Telegram.

    Parameters:
        bot_message : str
            Text message to send (used as the caption if a file is also
            sent, otherwise sent as a standalone text message).
        file_data : bytes, optional
            Raw bytes of a file (e.g. a JPG plot) to attach.
        file_name : str, optional
            Name to give the attached file (e.g. 'CMB_power_spectrum.jpg').
        caption : str, optional
            If given, overrides bot_message as the caption of the attached
            file. Kept separate for backwards compatibility.

    Returns:
        requests.Response
    """

    bot_token = '6271360875:AAGRCeeFGOpzw7EAnE_Koj8pwqLU1WldlAk'
    bot_chatID = '5908125195'

    if file_data is not None and file_name is not None: #send a picture
        send_file = 'https://api.telegram.org/bot' + bot_token + '/sendDocument?chat_id=' + bot_chatID
        file_bytes = io.BytesIO(file_data)
        data = {}
        text_caption = caption if caption is not None else bot_message
        if text_caption:
            data['caption'] = text_caption
        response = requests.post(send_file, files={'document': (file_name, file_bytes)}, data=data)
    else: #send a text message
        send_text = 'https://api.telegram.org/bot' + bot_token + '/sendMessage?chat_id=' + bot_chatID + '&parse_mode=Markdown&text=' + bot_message
        response = requests.get(send_text)

    return response



"""
bot_texter('The simulation has finished:') #for text

#For images
buffer = io.BytesIO()
plt.savefig(buffer, format='jpg')
buffer.seek(0)

image_base64 = base64.b64encode(buffer.getvalue()).decode()

bot_texter('Look', file_data=buffer.getvalue(), file_name='plots.jpg')
"""
