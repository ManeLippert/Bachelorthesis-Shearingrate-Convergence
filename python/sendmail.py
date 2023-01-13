import subprocess

def send_mail(recipient, subject, filename = None, body = None):

    recipient = recipient.encode("utf_8")
    
    subject = subject.replace(" ", "_")
    subject = subject.encode("utf_8")
    
    if body == None:
        body = "For futher information open attachment"
    
    body = body.encode("utf_8")
	
    if filename == None:
        process = subprocess.Popen(["ssh", "master", "/usr/bin/mailx", "-s", subject, recipient], stdin=subprocess.PIPE)
        process.communicate(body)

    else:
        attachmentPath = filename
        attachment = attachmentPath.encode("utf_8")
        
        process = subprocess.Popen(["ssh", "master", "/usr/bin/mailx", "-s", subject, "-a",  attachment, recipient], stdin=subprocess.PIPE)
        process.communicate(body)